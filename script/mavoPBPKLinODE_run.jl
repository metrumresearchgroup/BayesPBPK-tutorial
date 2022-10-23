#=
Version: October 18, 2022

Online supplement to the tutorial: 

Bayesian PBPK Modeling using R/Stan/Torsten and Julia/SciML/Turing.jl

This script runs the Bayesian PBPK linear ODE modeling in Julia/SciML/Turing.jl
=#

# activate pkg environment and load libraries
# libraries are referenced in the Project.toml file
cd(@__DIR__)
cd("..")
using Pkg; Pkg.activate(".") 
using CSV, DataFramesMeta, Chain, Random
using OrdinaryDiffEq, DiffEqCallbacks, Turing, Distributions, CategoricalArrays
using GlobalSensitivity, QuadGK
using Plots, StatsPlots, MCMCChains
using ExponentialUtilities
using LinearAlgebra
using Gadfly
import Cairo, Fontconfig

# set seed for reproducibility
Random.seed!(123)

# set paths
modDir = "model"
modName = "mavoPBPKLinODE"
tabDir = joinpath("deliv","table")
figDir = joinpath("deliv","figure")
modPath = mkpath(joinpath(modDir, string(modName, "_jl")))
figPath = mkpath(joinpath(figDir, string(modName, "_jl")))
tabPath = mkpath(joinpath(tabDir, string(modName, "_jl")))

# read mavoglurant Study A2121 data; we'll only use the first 20 subjects
dat_orig = CSV.read("data/Mavoglurant_A2121_nmpk.csv", DataFrame)
dat = dat_orig[dat_orig.ID .<= 812,:]
dat_obs = dat[dat.MDV .== 0,:]  # grab first 20 subjects ; remove missing obs

# load model
include(joinpath("..", modDir, string(modName, ".jl")));

# set simulation conditions
nSubject = length(unique(dat.ID))
doses = dat.AMT[dat.EVID .== 1,:]
rates = dat.RATE[dat.EVID .== 1,:]
durs = doses ./ rates
ncmt = 16
#u0 = zeros(ncmt)
tspan = (0.0,maximum(dat.TIME))
times = []
for i in unique(dat_obs.ID); push!(times, dat_obs.TIME[dat_obs.ID .== i]); end

# fixed parameters
wts = dat.WT[dat.EVID .== 1,:]
VVBs = (5.62 .* wts ./ 100) ./ 1.040  # volume of venous blood
BP = 0.61  # blood:plasma concentration ratio

# callback function for infusion stopping via zeroing k₀ parameter at end of infusion
function affect!(integrator)
    integrator.p[8] = integrator.p[8] == 0.0
end
cbs = []
for i in 1:length(unique(dat.ID)); push!(cbs, PresetTimeCallback([durs[i]], affect!)); end

# params to be estimated
CLint = exp(7.1)
KbBR = exp(1.1);
KbMU = exp(0.3);
KbAD = exp(2);
KbBO = exp(0.03);
KbRB = exp(0.3);
p = [CLint, KbBR, KbMU, KbAD, KbBO, KbRB, wts[1]]

# define the ODE problem
#prob = ODEProblem(PBPKODE!, u0, tspan, p)

#######

## Bayesian inference ##

## define Bayesian model function
@model function fitPBPK(data, nSubject, rates, times, wts, VVBs, BP) # data should be a Vector
    # priors
    ## informative priors for fixed effects and semi-informative for random effects
    σ ~ truncated(Cauchy(0, 0.5), 0.0, Inf)  # residual error
    ĈLint ~ LogNormal(7.1,0.25)  # intrinsic clearance
    KbBR ~ LogNormal(1.1,0.25)  # brain to plasma partition coefficient
    KbMU ~ LogNormal(0.3,0.25)  # muscle to plasma partition coefficient
    KbAD ~ LogNormal(2.0,0.25)  # adipose to plasma partition coefficient
    KbBO ~ LogNormal(0.03, 0.25)  # bone to plasma partition coefficient
    KbRB ~ LogNormal(0.3, 0.25)  # rest of body to plasma partition coefficient
    ω ~ truncated(Cauchy(0, 0.5), 0.0, 1.0)  # intersubject variability SD
    ηᵢ ~ filldist(Normal(0.0, 1.0), nSubject)  # individual ηs for random effects

    # individual params
    CLintᵢ = ĈLint .* exp.(ω .* ηᵢ)  # non-centered parameterization

    # initialize compartments with same type as params
    u0_tmp = eltype(KbBR).(zeros(ncmt))
    u0_tmp[15] = doses[1]
    rate_tmp = eltype(KbBR).(zeros(ncmt))
    rate_tmp[15] = rates[1]
    p1 = [CLintᵢ[1], KbBR, KbMU, KbAD, KbBO, KbRB, wts[1]]
    K_1 = PBPK(p1)
    K1 = DiffEqArrayOperator(K_1)
    time_reset = times[1] .- durs[1]
    prob1 = ODEProblem(K1, u0_tmp, (time_reset[1], time_reset[end]))
    
    # get preds for first two records for ID 1 that were measured before EOI
    u12 = eltype(K1).(zeros(ncmt, 2))
    u12 = K_1\ExponentialUtilities.expv_timestep(times[1][1:2], K1, rate_tmp) .- K1\rate_tmp

    # function to update ODE problem with newly sampled params
    u0_eoi = eltype(K1).(zeros(ncmt))
    function prob_func(prob,i,repeat)
        rate_tmp[15] = rates[i]
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts[i]]
        K_ = PBPK(ps)
        K = DiffEqArrayOperator(K_)
        u0_eoi = K\ExponentialUtilities.expv(durs[i], K, rate_tmp) - K\rate_tmp
        time_reset = times[i] .- durs[i]
        tspan_reset = (time_reset[1], time_reset[end])
        prob_tmp = ODEProblem(K, u0_eoi, tspan_reset, tstops=time_reset)
        return(prob_tmp)
    end

    # define an ensemble problem and simulate the population
    tmp_ensemble_prob = EnsembleProblem(prob1, prob_func=prob_func)
    # solving can be parallelized by replacing EnsembleSerial() with EnsembleThreads() but there is no need since chains will be parallelized later
    # and having nested parallelization might lead to non reproducible results
    tmp_ensemble_sol = solve(tmp_ensemble_prob, LinearExponential(), EnsembleSerial(), trajectories=nSubject) 

    # calculate predicted plasma concentration for the population and append in a vector 
    predicted = []
    for i in 1:nSubject
        times_tmp = times[i] .- durs[i]
        idx = findall(x -> x in times_tmp, tmp_ensemble_sol[i].t)
        tmp_sol = Array(tmp_ensemble_sol[i])[15,idx] ./ (VVBs[i]*BP/1000.0)
        append!(predicted, tmp_sol)
    end

    predicted[1:2] = u12[15,:] ./ (VVBs[1]*BP/1000.0)

    # likelihood
    for i = 1:length(predicted)
        data[i] ~ LogNormal(log.(max(predicted[i], 1e-12)), σ)
    end
end

# define the Bayesian model object
mod = fitPBPK(dat_obs.DV, nSubject, rates, times, wts, VVBs, BP)

# sample ; caution: sampling might take several hours to finish
# you can skip sampling and go directly to line 191 to load saved samples
## sampling conditions
nsampl = 250
nchains = 4
adapt_delta = .8

## sampling in serial 
#@time mcmcchains = sample(mod, NUTS(nsampl,adapt_delta), MCMCSerial(), nsampl, nchains)
#@time mcmcchains_prior = sample(mod, Prior(), MCMCSerial(), nsampl, nchains)  # parallel

## sampling with multithreading ; caution: results might not be reproducible 
@time mcmcchains = sample(mod, NUTS(nsampl,adapt_delta), MCMCThreads(), nsampl, nchains)
@time mcmcchains_prior = sample(mod, Prior(), MCMCThreads(), nsampl, nchains)  # parallel

## save mcmcchains
write(joinpath(modPath, string(modName, "chains.jls")), mcmcchains)
write(joinpath(modPath, string(modName, "chains_prior.jls")), mcmcchains_prior)

##uncomment to load saved chains
#mcmcchains = read(joinpath(modPath, string(modName, "chains.jls")), Chains)


#---# diagnostics #---#
# tables
summ, quant = describe(mcmcchains)
#summ, quant = describe(mcmcchains; q = [0.05, 0.25, 0.5, 0.75, 0.95])  # control quantiles in output

## summary
df_summ = DataFrame(summ)
CSV.write(joinpath(tabPath, "Summary.csv"), df_summ)

## quantiles
df_quant = DataFrame(quant)
CSV.write(joinpath(tabPath, "Quantiles.csv"), df_quant)

# plots
## trace plots ; we will split plots and join later for better view
#plot_chains = StatsPlots.plot(mcmcchains[:,1:8,:])  # mcmcchains[samples, params, chains]
plot_chains1 = StatsPlots.plot(mcmcchains[:,1:4,:])
plot_chains2 = StatsPlots.plot(mcmcchains[:,5:8,:])
plot_chains = Plots.plot(plot_chains1, plot_chains2, layout = (1,2))
savefig(plot_chains, joinpath(figPath, "MCMCTrace.pdf"))

## density plots
p_post = zeros(nsampl, 8, 4)
p_prior = deepcopy(p_post)

for i in 1:4; p_post[:,:,i] = Array(mcmcchains[:,1:8,1]); end
for i in 1:4; p_prior[:,:,i] = Array(mcmcchains_prior[:,1:8,1]); end

p_post_mean = mean(p_post, dims=3)[:,:,1]
p_prior_mean = mean(p_prior, dims=3)[:,:,1]

pars = summ[1:8,1]
dens_plots = []
for i in 1:8; p = density(p_post_mean[:,i], title=pars[i], label="Posterior"); density!(p_prior_mean[:,i], label="Prior"); push!(dens_plots, p); end

dens_plots[1] = Plots.plot(dens_plots[1], xlims=(0.0,0.6))

plot_dens = Plots.plot(dens_plots[1],
                        dens_plots[2],
                        dens_plots[3],
                        dens_plots[4],
                        dens_plots[5],
                        dens_plots[6],
                        dens_plots[7],
                        dens_plots[8], 
                        layout=grid(4,2),
                        size = (650,650))
savefig(plot_dens, joinpath(figPath, "DensPlots.pdf"))

#---# predictive checks #---#

#=
Predicitive checks can be conditioned on sampled parameters or on samples from a new population.
Here, we will do both.
=#

#--# conditional on chains #--#

# we will first create a vector of "missing" values to pass to the Bayesian fit function
# Turing will understand that the created model object is meant to simulate rather than fit
dat_missing = Vector{Missing}(missing, length(dat_obs.DV)) # vector of `missing`
mod_pred = fitPBPK(dat_missing, nSubject, rates, times, wts, VVBs, BP)
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

# data assembly to summarise and bring observed and predicted data together
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

## new etas for a new population
ηs = reshape(rand(Normal(0.0, 1.0), nSubject*nrow(df_params)), nrow(df_params), nSubject)

# run simulations for the new population
array_pred = Array{Float64}(undef, nrow(df_params), nrow(dat_obs))

for j in 1:nrow(df_params)
    KbBR = df_params[j,:KbBR]
    KbMU = df_params[j,:KbMU]
    KbAD = df_params[j,:KbAD]
    KbBO = df_params[j,:KbBO]
    KbRB = df_params[j,:KbRB]

    CLintᵢ = df_params[j,:ĈLint] .* exp.(df_params[j,:ω] .* ηs[j,:])

    # initialize compartments with same type as params
    u0_tmp = eltype(KbBR).(zeros(ncmt))
    u0_tmp[15] = doses[1]
    rate_tmp = eltype(KbBR).(zeros(ncmt))
    rate_tmp[15] = rates[1]
    p1 = [CLintᵢ[1], KbBR, KbMU, KbAD, KbBO, KbRB, wts[1]]
    K_1 = PBPK(p1)
    K1 = DiffEqArrayOperator(K_1)
    time_reset = times[1] .- durs[1]
    prob1 = ODEProblem(K1, u0_tmp, (time_reset[1], time_reset[end]))
    
    # get preds for first two records for ID 1 that were measured before EOI
    u12 = eltype(K1).(zeros(ncmt, 2))
    u12 = K_1\ExponentialUtilities.expv_timestep(times[1][1:2], K1, rate_tmp) .- K1\rate_tmp

    # function to update ODE problem with newly sampled params
    u0_eoi = eltype(K1).(zeros(ncmt))
    function prob_func(prob,i,repeat)
        rate_tmp[15] = rates[i]
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts[i]]
        K_ = PBPK(ps)
        K = DiffEqArrayOperator(K_)
        u0_eoi = K\ExponentialUtilities.expv(durs[i], K, rate_tmp) .- K\rate_tmp
        time_reset = times[i] .- durs[i]
        tspan_reset = (time_reset[1], time_reset[end])
        prob_tmp = ODEProblem(K, u0_eoi, tspan_reset, tstops=time_reset)
        return(prob_tmp)
    end

    # define an ensemble problem and simulate the population
    tmp_ensemble_prob = EnsembleProblem(prob1, prob_func=prob_func)
    tmp_ensemble_sol = solve(tmp_ensemble_prob, LinearExponential(), EnsembleSerial(), trajectories=nSubject) 

    predicted = []
    for i in 1:nSubject
        times_tmp = times[i] .- durs[i]
        idx = findall(x -> x in times_tmp, tmp_ensemble_sol[i].t)
        tmp_sol = Array(tmp_ensemble_sol[i])[15,idx] ./ (VVBs[i]*BP/1000.0)
        append!(predicted, tmp_sol)
    end

    predicted[1:2] = u12[15,:] ./ (VVBs[1]*BP/1000.0)

    array_pred[j, :] = rand.(LogNormal.(log.(max.(predicted, 1e-12)), df_params[j,:σ]))
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

# create the new population predictions dataframe
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

# conditioned on individual data
df_cCond = @chain begin
    df_vpc_pred
    groupby([:ID,:TIME])
    @transform(:loCond = quantile(:value, 0.05),
               :medCond = quantile(:value, 0.5),
               :hiCond = quantile(:value, 0.95))
    @select(:ID, :TIME, :DOSE, :loCond, :medCond, :hiCond)
    unique()
end

# conditioned on individuals in a new population
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





