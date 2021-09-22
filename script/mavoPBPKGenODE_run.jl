cd(@__DIR__)
using Pkg; Pkg.activate(".")
using CSV, DataFrames, Random
using OrdinaryDiffEq, DiffEqCallbacks, Turing, Distributions
using Plots, StatsPlots, MCMCChains

Random.seed!(123)

# read data
dat = CSV.read("data/Mavoglurant_A2121_nmpk.csv", DataFrame)
dat = dat[dat.ID .<= 812,:]
dat_obs = dat[dat.MDV .== 0,:]  # grab first 20 subjects ; remove missing obs

# load model
include("model/mavoPBPKGenODE.jl")

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
    ω ~ truncated(Cauchy(0, 0.5), 0.0, 2.0)
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
@time chains = mapreduce(c -> sample(mod, NUTS(250,.8), 250), chainscat, 1:4)  # serial
@time chains_prior = mapreduce(c -> sample(mod, Prior(), 250), chainscat, 1:4)  # serial

## multithreading
@time chains = sample(mod, NUTS(250,.8), MCMCThreads(), 250, 4)  # parallel
@time chains_prior = sample(mod, Prior(), MCMCThreads(), 250, 4)  # parallel

## save chains
save("model/mavoPBPKGenODEchains.jld","chains",chains)

## get results

## diagnostics
plot_chains = StatsPlots.plot(chains)  # chains[samples, params, chains]
#savefig(plot_chains, "BayesPopChains.pdf")

## predictive checks
dat_missing = Vector{Missing}(missing, length(dat_obs.DV)) # vector of `missing`
mod_pred = fitPBPK(dat_missing, prob)
pred = predict(mod_pred, chains)  # posterior
pred_prior = predict(mod_pred, chains_prior)

### summaries
summarystats(pred)
summ_pred, quant_pred = describe(pred)
summarystats(pred_prior)
summ_pred_prior, quant_pred_prior = describe(pred_prior)

#### save
#CSV.write("BayesPopSumm.csv", summ)
#CSV.write("BayesPopQuant.csv", quant)

### plot
plot_posteriorcheck = Gadfly.plot(x=data_obs.TIME, y=data_obs.DV, Geom.point, Theme(background_color = "white"), Guide.xlabel("Time"), Guide.ylabel("Concentration"), Guide.title("Posterior predictive check"),
    layer(x=data_obs.TIME, y=quant_pred[:,4], Geom.line),
    layer(x=data_obs.TIME, ymin=quant_pred[:,2], ymax=quant_pred[:,6], Geom.ribbon))

plot_priorcheck = Gadfly.plot(x=data_obs.TIME, y=data_obs.DV, Geom.point, Theme(background_color = "white"), Guide.xlabel("Time"), Guide.ylabel("Concentration"), Guide.title("Prior predictive check"),
    layer(x=data_obs.TIME, y=quant_pred_prior[:,4], Geom.line),
    layer(x=data_obs.TIME, ymin=quant_pred_prior[:,2], ymax=quant_pred_prior[:,6], Geom.ribbon))

hstack(plot_priorcheck, plot_posteriorcheck)

