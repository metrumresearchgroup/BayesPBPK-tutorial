cd(@__DIR__)
using Pkg; Pkg.activate(".")
using CSV, DataFramesMeta, Chain, Random
using OrdinaryDiffEq, DiffEqCallbacks, Turing, Distributions, CategoricalArrays
using GlobalSensitivity, QuadGK
using Plots, StatsPlots, MCMCChains
using Gadfly
import Cairo, Fontconfig

Random.seed!(1234)

# paths
modDir = "model"
modName = "mavoPBPKGenODE"
tabDir = joinpath("deliv","table")
figDir = joinpath("deliv","figure")
figPath = mkpath(joinpath(figDir, string(modName, "_jl")))
tabPath = mkpath(joinpath(tabDir, string(modName, "_jl")))

# read data
dat_orig = CSV.read("data/Mavoglurant_A2121_nmpk.csv", DataFrame)
dat = dat_orig[dat_orig.ID .<= 812,:]
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

########

## sensitivity analysis ##
## create function that takes in paramers and returns endpoints for sensitivity
p_sens = p[1:6]
f_globsens = function(p_sens)
    p_all = [p_sens;[wts[1],rates[1]]]
    tmp_prob = remake(prob, p = p_all)
    tmp_sol = solve(tmp_prob, Tsit5(), save_idxs=[15], callback=cbs[1])
    auc, err = quadgk(tmp_sol, 0.0, 48.0)
    return(auc[1])
end

### bounds
bounds = [[1000.0,1500.0],[1.0,10.0],[1.0,10.0],[1.0,10.0],[1.0,10.0],[1.0,10.0]]

#### Sobol
@time s = GlobalSensitivity.gsa(f_globsens, Sobol(), bounds, N=1000)

plot_sens_total = Plots.bar(["CLint","KbBR","KbMU","KbAD","KbBO","KbRB"], s.ST, title="Total Order Indices", ylabel="Index", legend=false)
Plots.hline!([0.05], linestyle=:dash)
plot_sens_single = Plots.bar(["CLint","KbBR","KbMU","KbAD","KbBO","KbRB"], s.S1, title="First Order Indices", ylabel="Index", legend=false)
Plots.hline!([0.05], linestyle=:dash)

plot_sens = Plots.plot(plot_sens_single, plot_sens_total, xrotation = 45)
savefig(plot_sens, joinpath(figPath, "sensitivity.pdf"))

#######

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
#@time mcmcchains = mapreduce(c -> sample(mod, NUTS(250,.8), 250), chainscat, 1:4)  # serial
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
#summ, quant = describe(mcmcchains; q = [0.05, 0.25, 0.5, 0.75, 0.95])

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

#--# conditional on chains #--#

dat_missing = Vector{Missing}(missing, length(dat_obs.DV)) # vector of `missing`
mod_pred = fitPBPK(dat_missing, prob, nSubject, rates, times, wts, cbs, VVBs, BP)
#pred = predict(mod_pred, mcmcchains)  # posterior ; conditioned on each sample in chains
pred = predict(mod_pred, mcmcchains, include_all=false)  # include_all = false means sampling new !!
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
bins = [0, 1, 2, 3, 4, 6, 8, 10, 20, 30, 40, 50]
labels = string.(1:length(bins) - 1)

## observed
df_vpc_obs = @chain begin
    dat_obs
    @select(:ID, :TIME, :DV, :DOSE)
    @transform(:DNDV = :DV ./ :DOSE,
               :bins = cut(:TIME, bins, labels = labels))
    
    groupby(:bins)
    @transform(:lo = quantile(:DNDV, 0.05),
               :med = quantile(:DNDV, 0.5),
               :hi = quantile(:DNDV, 0.95))
end
df_vpc_obs2 = @orderby(unique(@select(df_vpc_obs, :TIME, :bins, :lo, :med, :hi)), :TIME)

## predicted
df_pred = DataFrame(pred)
df_vpc_pred = @chain begin
    df_pred
    DataFramesMeta.stack(3:ncol(df_pred))
    @orderby(:iteration, :chain)
    hcat(select(repeat(dat_obs, 1000), [:ID,:TIME,:DOSE]))
    @transform(:DNDV = :value ./ :DOSE,
    :bins = cut(:TIME, bins, labels = labels))

    groupby([:iteration, :chain, :bins])
    @transform(:lo = quantile(:DNDV, 0.05),
               :med = quantile(:DNDV, 0.5),
               :hi = quantile(:DNDV, 0.95))
    
    groupby(:bins)
    @transform(:loLo = quantile(:lo, 0.025),
               :medLo = quantile(:lo, 0.5),
               :hiLo = quantile(:lo, 0.975),
               :loMed = quantile(:med, 0.025),
               :medMed = quantile(:med, 0.5),
               :hiMed = quantile(:med, 0.975),
               :loHi = quantile(:hi, 0.025),
               :medHi = quantile(:hi, 0.5),
               :hiHi = quantile(:hi, 0.975))
end

df_vpc_pred2 = @orderby(unique(df_vpc_pred[!,[6;13:21]]), :TIME)

### plot
dat_obs2 = @transform(dat_obs, :DNDV = :DV ./ :DOSE)

set_default_plot_size(17cm, 12cm)

plot_ppc = Gadfly.plot(x=dat_obs2.TIME, y=dat_obs2.DNDV, Geom.point, Scale.y_log10, Theme(background_color="white", default_color="black"), alpha=[0.2], Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant dose-normalized concentration (ng/mL/mg)", orientation=:vertical),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.med, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.lo, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.hi, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_pred2.TIME, ymin=df_vpc_pred2.loMed, ymax=df_vpc_pred2.hiMed, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.8]),
    layer(x=df_vpc_pred2.TIME, ymin=df_vpc_pred2.loLo, ymax=df_vpc_pred2.hiLo, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]),
    layer(x=df_vpc_pred2.TIME, ymin=df_vpc_pred2.loHi, ymax=df_vpc_pred2.hiHi, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]))

plot_tmp = PDF(joinpath(figPath, "PPCCond.pdf"), 17cm, 12cm)
draw(plot_tmp, plot_ppc)


#--# new population #--#

df_params = DataFrame(mcmcchains)[:,3:10]

## new etas
ηs = reshape(rand(Normal(0.0, 1.0), nSubject*nrow(df_params1)), nrow(df_params1), nSubject)

array_pred = Array{Float64}(undef, nrow(df_params), nrow(dat_obs))

for j in 1:nrow(df_params)
    KbBR = df_params[j,:KbBR]
    KbMU = df_params[j,:KbMU]
    KbAD = df_params[j,:KbAD]
    KbBO = df_params[j,:KbBO]
    KbRB = df_params[j,:KbRB]

    CLintᵢ = df_params[j,:ĈLint] .* exp.(df_params[j,:ω] .* ηs[j,:])

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

    #array_pred[j, :] = predicted
    array_pred[j, :] = rand.(LogNormal.(log.(predicted), df_params[j,:σ]))
end

df_pred_new = DataFrame(array_pred, :auto)
@transform!(df_pred_new, :iteration = 1:size(array_pred)[1])

df_vpc_pred_new = @chain begin
    df_pred_new
    DataFramesMeta.stack(1:268)
    @orderby(:iteration)
    hcat(select(repeat(dat_obs, 1000), [:ID,:TIME,:DOSE]))
    @transform(:DNDV = :value ./ :DOSE,
    :bins = cut(:TIME, bins, labels = labels))

    groupby([:iteration, :bins])
    @transform(:lo = quantile(:DNDV, 0.05),
               :med = quantile(:DNDV, 0.5),
               :hi = quantile(:DNDV, 0.95))
    
    groupby(:bins)
    @transform(:loLo = quantile(:lo, 0.025),
               :medLo = quantile(:lo, 0.5),
               :hiLo = quantile(:lo, 0.975),
               :loMed = quantile(:med, 0.025),
               :medMed = quantile(:med, 0.5),
               :hiMed = quantile(:med, 0.975),
               :loHi = quantile(:hi, 0.025),
               :medHi = quantile(:hi, 0.5),
               :hiHi = quantile(:hi, 0.975))
end

df_vpc_pred_new2 = @orderby(unique(df_vpc_pred_new[!,[5;12:20]]), :TIME)

### plot
#dat_obs2 = @transform(dat_obs, :DNDV = :DV ./ :DOSE)

set_default_plot_size(17cm, 12cm)

plot_ppc_new = Gadfly.plot(x=dat_obs2.TIME, y=dat_obs2.DNDV, Geom.point, Scale.y_log10, Theme(background_color="white", default_color="black"), alpha=[0.2], Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant dose-normalized concentration (ng/mL/mg)", orientation=:vertical),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.med, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.lo, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_obs.TIME, y=df_vpc_obs.hi, Geom.line, Theme(default_color="black")),
    layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loMed, ymax=df_vpc_pred_new2.hiMed, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.8]),
    layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loLo, ymax=df_vpc_pred_new2.hiLo, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]),
    layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loHi, ymax=df_vpc_pred_new2.hiHi, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]))

plot_tmp = PDF(joinpath(figPath, "PPCPred.pdf"), 17cm, 12cm)
draw(plot_tmp, plot_ppc_new)

#######

## simulation ##

# get wts for 500 individuals
wts_sim = @chain begin 
    dat_orig 
    @subset(:EVID .== 1) 
    unique(:ID)
    @select(:WT)
    Array
end
wts_sim = sample(wts_sim, 500, replace=true)
vvbs_sim = (5.62 .* wts_sim ./ 100) ./ 1.040

# get 500 replicates of inferred parameters
df_params1 = DataFrame(mcmcchains)[:,3:10]
df_params2 = df_params1[sample(1:nrow(df_params1), 500, replace=false),:]

# simulate a 50 mg single and multiple daily doses



#=
df_params2 = @chain begin
    df_params1
    @transform(:ηᵢ = rand(Normal(), nrow(df_params1)),
               :WT = rand(wts, nrow(df_params1)))
    @transform(:CLintᵢ = :ĈLint .* exp.(:ω .* :ηᵢ),
               :VVB = (5.62 .* :WT ./ 100) ./ 1.040,
               :RATE = 50.0 / (10.0/60.0))  # 50 mg infused over 10 min
end

# run simulation
## infusion callback
cb = PresetTimeCallback([10.0/60.0], affect!)

### define new problem for simulation
function prob_func_sim(prob,i,repeat)
    ps = [df_params2.CLintᵢ[i], df_params2.KbBR[i], df_params2.KbMU[i], df_params2.KbAD[i], df_params2.KbBO[i], df_params2.KbRB[i], df_params2.WT[i], df_params2.RATE[i]]
    remake(prob, p=ps, tspan=(0.0,48.0))
end

### simulate
sim_ensemble_prob = EnsembleProblem(prob, prob_func=prob_func_sim)
sim_ensemble_sol = solve(sim_ensemble_prob, Tsit5(), save_idxs=[15], callback=cb, saveat=0.1, trajectories=nrow(df_params2)) 

# calc stats
df_sim1 = []
for i in 1:nrow(df_params1)
    df_tmp = @chain begin
        df_tmp = DataFrame(sim_ensemble_sol[i])
        rename(:timestamp => :time, :value1 => :cent)
        @transform!(:ID = i,
                    :cent = :cent / (df_params2.VVB[i]*BP/1000.0)) 
    end
    push!(df_sim1, df_tmp)
end

df_sim2 = @chain begin
    vcat(df_sim1...)
    groupby(:time)
    @transform(:med = quantile(:cent, 0.5),
               :lo = quantile(:cent, 0.05),
               :hi = quantile(:cent, 0.95))
    @subset(:ID .== 1)
    unique
end

set_default_plot_size(17cm, 12cm)

plot_sim = Gadfly.plot(x=df_sim2.time, y=df_sim2.med, Geom.line, Scale.y_log10, Theme(background_color="white", default_color="black"), Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant concentration (ng/mL)", orientation=:vertical),
    layer(x=df_sim2.time, ymin=df_sim2.lo, ymax=df_sim2.hi, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]))

plot_tmp = PDF(joinpath(figPath, "sim.pdf"), 17cm, 12cm)
draw(plot_tmp, plot_sim)
=#


