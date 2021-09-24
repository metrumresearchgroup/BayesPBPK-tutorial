cd(@__DIR__)
using Pkg; Pkg.activate(".")
using CSV, DataFrames, Random
using OrdinaryDiffEq, DiffEqCallbacks, Turing, Distributions
using Plots, StatsPlots, MCMCChains
using Gadfly

Random.seed!(1234)

# paths
modDir = "model"
modName = "mavoPBPKGenODE"
tabDir = joinpath("deliv","table")
figDir = joinpath("deliv","figure")
figPath = mkpath(joinpath(figDir, string(modName, "_jl")))
tabPath = mkpath(joinpath(tabDir, string(modName, "_jl")))

# read data
dat = CSV.read("data/Mavoglurant_A2121_nmpk.csv", DataFrame)
dat = dat[dat.ID .<= 812,:]
dat_obs = dat[dat.MDV .== 0,:]  # grab first 20 subjects ; remove missing obs

# load model
include(joinpath(modDir, string(modName, ".jl")))

# conditions
nSubject = length(unique(dat.ID))
doses = dat.AMT[dat.EVID .== 1,:]
rates = dat.RATE[dat.EVID .== 1,:]
durs = doses ./ rates

ncmt = 16
u0 = zeros(ncmt)
tspan = (0.0,maximum(dat.TIME))
times = []
for i in unique(dat_obs.ID); push!(times, dat_obs.TIME[dat_obs.ID .== i]); end

# fixed parameters
wts = dat.WT[dat.EVID .== 1,:]
VVBs = (5.62 .* wts ./ 100) ./ 1.040
BP = 0.61

# callback function for infusion stopping
function affect!(integrator)
    integrator.p[8] = integrator.p[8] == 0.0
end
cbs = []
for i in 1:length(unique(dat.ID)); push!(cbs, PresetTimeCallback([durs[i]], affect!)); end

# variable params
CLint = exp(7.1)
KbBR = exp(1.1);
KbMU = exp(0.3);
KbAD = exp(2);
KbBO = exp(0.03);
KbRB = exp(0.3);
p = [CLint, KbBR, KbMU, KbAD, KbBO, KbRB, wts[1], rates[1]]

# define problem
prob = ODEProblem(PBPKODE!, u0, tspan, p)

## Bayesian inference ##
@model function fitPBPK(data, prob, nSubject, rates, times, wts, cbs, VVBs, BP) # data should be a Vector
    # priors
    σ ~ truncated(Cauchy(0, 0.5), 0.0, Inf) # residual error
    ĈLint ~ LogNormal(7.1,0.25)
    KbBR ~ LogNormal(1.1,0.25)
    KbMU ~ LogNormal(0.3,0.25)
    KbAD ~ LogNormal(2.0,0.25)
    KbBO ~ LogNormal(0.03, 0.25)
    KbRB ~ LogNormal(0.3, 0.25)
    ω ~ truncated(Cauchy(0, 0.5), 0.0, 1.0)
    ηᵢ ~ filldist(Normal(0.0, 1.0), nSubject)

    # individual params
    CLintᵢ = ĈLint .* exp.(ω .* ηᵢ)  # non-centered parameterization

    # simulate
    function prob_func(prob,i,repeat)
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts[i], rates[i]]
        remake(prob, p=ps, saveat=times[i], callback=cbs[i])
    end

    tmp_ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    tmp_ensemble_sol = solve(tmp_ensemble_prob, Tsit5(), trajectories=nSubject) 

    predicted = []
    for i in 1:nSubject
        times_tmp = times[i]
        idx = findall(x -> x in times_tmp, tmp_ensemble_sol[i].t)
        tmp_sol = Array(tmp_ensemble_sol[i])[15,idx] ./ (VVBs[i]*BP/1000.0)
        append!(predicted, tmp_sol)
    end

    # likelihood
    for i = 1:length(predicted)
        data[i] ~ LogNormal(log.(predicted[i]), σ)
    end
end

mod = fitPBPK(dat_obs.DV, prob, nSubject, rates, times, wts, cbs, VVBs, BP)

# run 
## serial 
@time mcmcchains = mapreduce(c -> sample(mod, NUTS(250,.8), 250), chainscat, 1:4)  # serial
#@time mcmcchains_prior = mapreduce(c -> sample(mod, Prior(), 250), chainscat, 1:4)  # serial

## multithreading
@time mcmcchains = sample(mod, NUTS(250,.8), MCMCThreads(), 250, 4)  # parallel
#@time mcmcchains_prior = sample(mod, Prior(), MCMCThreads(), 250, 4)  # parallel

## save mcmcchains
#write(joinpath(modDir, string(modName, "chains.jls")), mcmcchains)

##load saved chains
#mcmcchains = read(joinpath(modDir, string(modName, "chains.jls")), Chains)


#---# diagnostics #---#
# tables
summ, quant = describe(mcmcchains)

## summary
df_summ = DataFrame(summ)
CSV.write(joinpath(tabPath, "Summary.csv"), df_summ)

## quantiles
df_quant = DataFrame(quant)
CSV.write(joinpath(tabPath, "Quantiles.csv"), df_quant)

# plots
## trace plots
#plot_chains = StatsPlots.plot(mcmcchains[:,1:8,:])  # mcmcchains[samples, params, chains]
plot_chains1 = StatsPlots.plot(mcmcchains[:,1:4,:])
plot_chains2 = StatsPlots.plot(mcmcchains[:,5:8,:])
plot_chains = Plots.plot(plot_chains1, plot_chains2, layout = (1,2))
savefig(plot_chains, joinpath(figPath, "MCMCTrace.pdf"))

#---# predictive checks #---#
dat_missing = Vector{Missing}(missing, length(dat_obs.DV)) # vector of `missing`
mod_pred = fitPBPK(dat_missing, prob, nSubject, rates, times, wts, cbs, VVBs, BP)
pred = predict(mod_pred, mcmcchains)  # posterior
#pred_prior = predict(mod_pred, mcmcchains_prior)

### predictive checks summaries
#summarystats(pred)
summ_pred, quant_pred = describe(pred)
#summarystats(pred_prior)
#summ_pred_prior, quant_pred_prior = describe(pred_prior)

#### save
#CSV.write(joinpath(tabPath, "Summary_PPC.csv"), summ_pred)
#CSV.write(joinpath(tabPath, "Quantiles_PPC.csv"), quant_pred)

# data assembly
bins = [0, 2, 4, 6, 8, 10, 20, 30, 40, 50]
labels = string.(1:length(bins) - 1)

## observed
dat_obs_vpc1 = DataFrame(ID = dat_obs.ID, TIME = dat_obs.TIME, DV = dat_obs.DV, DOSE = dat_obs.DOSE)
dat_obs_vpc2 = transform(dat_obs_vpc1, [:DV, :DOSE] => ByRow((x,y) -> x/y) => :DNDV)
dat_obs_vpc2.bins = cut(dat_obs_vpc2.TIME, bins, labels=labels)

dat_obs_vpc3 = transform(groupby(dat_obs_vpc2, :bins), :DNDV => x -> quantile(x, 0.05))
rename!(dat_obs_vpc3, :DNDV_function => :lo)
dat_obs_vpc3 = transform(groupby(dat_obs_vpc3, :bins), :DNDV => x -> quantile(x, 0.5))
rename!(dat_obs_vpc3, :DNDV_function => :med)
dat_obs_vpc3 = transform(groupby(dat_obs_vpc3, :bins), :DNDV => x -> quantile(x, 0.95))
rename!(dat_obs_vpc3, :DNDV_function => :hi)

## predicted
df_pred1 = DataFrame(pred)
df_pred2 = DataFrames.stack(df_pred1, 3:270)
sort!(df_pred2, [:iteration,:chain])
dat_obs_rep = select(repeat(dat_obs, 1000), [:ID,:TIME,:DOSE])
df_pred3 = hcat(dat_obs_rep, df_pred2)
df_pred4 = transform(df_pred3, [:value, :DOSE] => ByRow((x,y) -> x/y) => :DNDV)
df_pred4.bins = cut(df_pred4.TIME, bins, labels=labels)

df_pred5 = transform(groupby(df_pred4, [:iteration,:bins]), :DNDV => x -> quantile(x, 0.05))
rename!(df_pred5, :DNDV_function => :lo)
df_pred5 = transform(groupby(df_pred5, [:iteration,:bins]), :DNDV => x -> quantile(x, 0.5))
rename!(df_pred5, :DNDV_function => :med)
df_pred5 = transform(groupby(df_pred5, [:iteration,:bins]), :DNDV => x -> quantile(x, 0.95))
rename!(df_pred5, :DNDV_function => :hi)

df_pred5 = transform(groupby(df_pred5, [:bins]), :lo => x -> quantile(x, 0.05))
rename!(df_pred5, :lo_function => :loLo)
df_pred5 = transform(groupby(df_pred5, [:bins]), :lo => x -> quantile(x, 0.5))
rename!(df_pred5, :lo_function => :medLo)
df_pred5 = transform(groupby(df_pred5, [:bins]), :lo => x -> quantile(x, 0.95))
rename!(df_pred5, :lo_function => :hiLo)

df_pred5 = transform(groupby(df_pred5, [:bins]), :med => x -> quantile(x, 0.05))
rename!(df_pred5, :med_function => :loMed)
df_pred5 = transform(groupby(df_pred5, [:bins]), :med => x -> quantile(x, 0.5))
rename!(df_pred5, :med_function => :medMed)
df_pred5 = transform(groupby(df_pred5, [:bins]), :med => x -> quantile(x, 0.95))
rename!(df_pred5, :med_function => :hiMed)

df_pred5 = transform(groupby(df_pred5, [:bins]), :hi => x -> quantile(x, 0.05))
rename!(df_pred5, :hi_function => :loHi)
df_pred5 = transform(groupby(df_pred5, [:bins]), :hi => x -> quantile(x, 0.5))
rename!(df_pred5, :hi_function => :medHi)
df_pred5 = transform(groupby(df_pred5, [:bins]), :hi => x -> quantile(x, 0.95))
rename!(df_pred5, :hi_function => :hiHi)

df_pred6 = subset(df_pred5, :iteration => ByRow(==(1)), :chain => ByRow(==(1)))

### plot
set_default_plot_size(17cm, 12cm)
plot_posteriorcheck = Gadfly.plot(x=dat_obs_vpc3.TIME, y=dat_obs_vpc3.DNDV, Geom.point, Scale.y_log10, Theme(background_color = "white"), Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant dose-normalized concentration (ng/mL/mg)", orientation=:vertical),
    layer(x=dat_obs_vpc3.TIME, y=dat_obs_vpc3.med, Geom.line),
    layer(x=dat_obs_vpc3.TIME, y=dat_obs_vpc3.lo, Geom.line),
    layer(x=dat_obs_vpc3.TIME, y=dat_obs_vpc3.hi, Geom.line),
    layer(x=df_pred6.TIME, ymin=df_pred6.loLo, ymax=df_pred6.hiLo, Geom.ribbon),
    layer(x=df_pred6.TIME, ymin=df_pred6.loMed, ymax=df_pred6.hiMed, Geom.ribbon),
    layer(x=df_pred6.TIME, ymin=df_pred6.loHi, ymax=df_pred6.hiHi, Geom.ribbon))
    
p = PDF(joinpath(figPath, "PPC.pdf"), 17cm, 12cm)
draw(p, plot_posteriorcheck)

