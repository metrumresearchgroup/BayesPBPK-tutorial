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
modPath = mkpath(joinpath(modDir, string(modName, "_jl")))
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

#=
## rhat
df_tmp = @orderby(DataFrame(summ)[1:8,:], :rhat)
plot_rhat = Plots.bar(string.(df_tmp.parameters), df_tmp.rhat, orientation=:h, legend = false, xlim=[0.99,1.05], xticks=[1.0,1.05], xlabel = "R̂")
Plots.vline!([1.0,1.05], linestyle=[:solid,:dash])
savefig(plot_rhat, joinpath(figPath, "rhat.pdf"))

## neff
df_tmp = @orderby(@transform!(df_tmp, :neff = :ess ./ 1000.0), :neff)
plot_neff = Plots.bar(string.(df_tmp.parameters), df_tmp.neff, orientation=:h, legend = false, xlim=[0.0,2.5], xticks=[0.0:0.25:2.5;], xlabel = "Neff/N")
Plots.vline!([0.1,0.5,1.0], linestyle=:dash)
savefig(plot_neff, joinpath(figPath, "neff.pdf"))
=#

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

# save CSV
#CSV.write(joinpath(modPath, "df_params.csv"), df_params)

## new etas
ηs = reshape(rand(Normal(0.0, 1.0), nSubject*nrow(df_params)), nrow(df_params), nSubject)

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

    array_pred[j, :] = rand.(LogNormal.(log.(predicted), df_params[j,:σ]))
end

df_pred_new = DataFrame(array_pred, :auto)
@transform!(df_pred_new, :iteration = 1:size(array_pred)[1])

# save version for R's vpc
df_pred_new2 = @chain begin
    df_pred_new
    DataFramesMeta.stack(1:268)
    @orderby(:iteration)
    hcat(select(repeat(dat_obs, 1000), [:ID,:TIME,:DOSE]))
    @transform(:DNDV = :value ./ :DOSE)
end

# save CSV
#CSV.write(joinpath(modPath, "df_pred.csv"), df_pred_new2)

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

#--# individual plots #--#

df_cObs = @chain begin
    df_vpc_obs
    @select(:ID,:TIME,:DV,:DOSE)
end

df_cCond = @chain begin
    df_vpc_pred
    groupby([:ID,:TIME])
    @transform(:loCond = quantile(:value, 0.05),
               :medCond = quantile(:value, 0.5),
               :hiCond = quantile(:value, 0.95))
    @select(:ID, :TIME, :DOSE, :loCond, :medCond, :hiCond)
    unique()
end

df_cPred = @chain begin
    df_vpc_pred_new
    groupby([:ID,:TIME])
    @transform(:loPred = quantile(:value, 0.05),
               :medPred = quantile(:value, 0.5),
               :hiPred = quantile(:value, 0.95))
    @select(:ID, :TIME, :DOSE, :loPred, :medPred, :hiPred)
    unique()
end

# join all data
df_cAll = hcat(df_cObs, 
               @select(df_cCond, :loCond, :medCond, :hiCond),
               @select(df_cPred, :loPred, :medPred, :hiPred))

# save CSV for plotting in R
CSV.write(joinpath(modPath, "df_ind.csv"), df_cAll)

#=
# plot
plot_ind = Gadfly.plot(df_call, x=:TIME, y=:DV, xgroup=:ID, Geom.subplot_grid(Geom.point),
                       Theme(background_color="white", default_color="black"), 
                       Scale.y_log10) 
plot_ind = Gadfly.plot(df_call, x=:TIME, Theme(background_color="white", default_color="black"), Scale.y_log10, 
    layer(y=:DV, Geom.point))


plot_ind = Gadfly.plot(x=dat_obs2.TIME, y=dat_obs2.DNDV, Geom.point, Scale.y_log10, Theme(background_color="white", default_color="black"), alpha=[0.2], Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant dose-normalized concentration (ng/mL/mg)", orientation=:vertical),
layer(x=df_vpc_obs.TIME, y=df_vpc_obs.med, Geom.line, Theme(default_color="black")),
layer(x=df_vpc_obs.TIME, y=df_vpc_obs.lo, Geom.line, Theme(default_color="black")),
layer(x=df_vpc_obs.TIME, y=df_vpc_obs.hi, Geom.line, Theme(default_color="black")),
layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loMed, ymax=df_vpc_pred_new2.hiMed, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.8]),
layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loLo, ymax=df_vpc_pred_new2.hiLo, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]),
layer(x=df_vpc_pred_new2.TIME, ymin=df_vpc_pred_new2.loHi, ymax=df_vpc_pred_new2.hiHi, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]))
=#

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
df_params_sim1 = DataFrame(mcmcchains)[:,3:10]
df_params_sim2 = df_params_sim1[sample(1:nrow(df_params_sim1), 500, replace=false),:]

## new etas
ηs_sim = reshape(rand(Normal(0.0, 1.0), 500*500), 500, 500)

# create array to hold results
times_sim = [0.0:0.1:48.0;]
array_pred_sim = Array{Float64}(undef, 500, length(times_sim)*500)

# simulate a 50 mg single
## infusion callback
dose = 50.0
cb = PresetTimeCallback([10.0/60.0], affect!)

for j in 1:nrow(df_params_sim2)
    KbBR = df_params_sim2[j,:KbBR]
    KbMU = df_params_sim2[j,:KbMU]
    KbAD = df_params_sim2[j,:KbAD]
    KbBO = df_params_sim2[j,:KbBO]
    KbRB = df_params_sim2[j,:KbRB]

    CLintᵢ = df_params_sim2[j,:ĈLint] .* exp.(df_params_sim2[j,:ω] .* ηs_sim[j,:])

    # simulate
    function prob_func_sim(prob,i,repeat)
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts_sim[i], 300.0]
        remake(prob, p=ps, tspan=(0.0,48.0))
    end
    
    sim_ensemble_prob = EnsembleProblem(prob, prob_func=prob_func_sim)
    sim_ensemble_sol = solve(sim_ensemble_prob, Tsit5(), save_idxs=[15], callback=cb, saveat=times_sim, trajectories=nrow(df_params_sim2))

    predicted = []
    for i in 1:500
        idx = findall(x -> x in times_sim, sim_ensemble_sol[i].t)
        tmp_sol = Array(sim_ensemble_sol[i])[1,idx] ./ (vvbs_sim[i]*BP/1000.0)
        append!(predicted, tmp_sol)
    end

    array_pred_sim[j, :] = rand.(LogNormal.(log.(predicted), df_params_sim2[j,:σ]))
end

# get stats

df_pred_sim = DataFrame(array_pred_sim, :auto)
@transform!(df_pred_sim, :iteration = 1:size(array_pred_sim)[1])

df_pred_sim2 = @chain begin
    df_pred_sim
    DataFramesMeta.stack(1:240500)
    @orderby(:iteration)
    @transform(:ID = repeat(repeat(1:500,inner=length(times_sim)), 500),
               :TIME = repeat(times_sim, 500*500),
               :DOSE = dose)
    @transform(:DNDV = :value ./ :DOSE)

    groupby([:iteration, :TIME])
    @transform(:lo = quantile(:value, 0.05),
               :med = quantile(:value, 0.5),
               :hi = quantile(:value, 0.95))
    
    groupby(:TIME)
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

df_pred_summ = @orderby(unique(df_pred_sim2[!,[5;11:19]]), :TIME)

# save CSV
#CSV.write(joinpath(modPath, "df_pred.csv"), df_pred_new2)

### plot
#dat_obs2 = @transform(dat_obs, :DNDV = :DV ./ :DOSE)

set_default_plot_size(17cm, 12cm)

plot_pred_summ = Gadfly.plot(x=df_pred_summ.TIME, ymin=df_pred_summ.loMed, ymax=df_pred_summ.hiMed, Geom.ribbon, Scale.y_log10, Theme(default_color="deepskyblue", background_color="white"), alpha=[0.8], Guide.xlabel("Time (h)"), Guide.ylabel("Mavoglurant concentration (ng/mL)", orientation=:vertical),
    layer(x=df_pred_summ.TIME, ymin=df_pred_summ.loLo, ymax=df_pred_summ.hiLo, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]),
    layer(x=df_pred_summ.TIME, ymin=df_pred_summ.loHi, ymax=df_pred_summ.hiHi, Geom.ribbon, Theme(default_color="deepskyblue"), alpha=[0.5]),
    layer(x=df_pred_summ.TIME, y=df_pred_summ.medMed, Geom.line, Theme(default_color="black")),
    layer(x=df_pred_summ.TIME, y=df_pred_summ.medLo, Geom.line, Theme(default_color="black")),
    layer(x=df_pred_summ.TIME, y=df_pred_summ.medHi, Geom.line, Theme(default_color="black")))

plot_tmp = PDF(joinpath(figPath, "SimPred.pdf"), 17cm, 12cm)
draw(plot_tmp, plot_pred_summ)



