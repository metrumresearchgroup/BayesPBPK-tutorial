
cd(@__DIR__)
cd("..")
using Pkg; Pkg.activate(".")
using CSV, DataFrames, Random
using OrdinaryDiffEq, DiffEqCallbacks, Turing, Distributions, ExponentialAction
using Plots, StatsPlots, MCMCChains
using Gadfly

Random.seed!(123456)

# paths
modDir = "model"
modName = "mavoPBPKLinODE"
tabDir = joinpath("deliv","table")
figDir = joinpath("deliv","figure")
modPath = mkpath(joinpath(modDir, string(modName, "_jl")))
figPath = mkpath(joinpath(figDir, string(modName, "_jl")))
tabPath = mkpath(joinpath(tabDir, string(modName, "_jl")))

# read data
dat = CSV.read("data/Mavoglurant_A2121_nmpk.csv", DataFrame)
dat = dat[dat.ID .<= 812,:]
dat_obs = dat[dat.MDV .== 0,:]  # grab first 20 subjects ; remove missing obs

# load model
include(joinpath("..", modDir, string(modName, ".jl")));

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

# housekeeping
ut = zeros(ncmt,nrow(dat_obs))
rate_tmp = zeros(ncmt)
lens = []
for i in 1:length(times); push!(lens, length(times[i])); end
sum_lens = cumsum(lens)
starts = [1;sum_lens[1:end - 1] .+ 1]
ends = deepcopy(sum_lens)

## Bayesian inference ##
@model function fitPBPK(data, nSubject, rates, times, durs, wts, VVBs, BP, starts, ends, ncmt) # data should be a Vector
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

    # solve
    tmp_rate = zeros(ncmt)
    tmp_rate[15] = rates[1]
    predicted = []

    ps_tmp = [CLintᵢ[1], KbBR, KbMU, KbAD, KbBO, KbRB, wts[1]]
    K_tmp = PBPK(ps_tmp)
    u0_eoi_tmp = inv(K_tmp)*ExponentialAction.expv(durs[1], K_tmp, tmp_rate) - inv(K_tmp)*tmp_rate
    ut = zeros(Base.promote_eltype(times[1], K_tmp, u0_eoi_tmp), length(u0_eoi_tmp), length(data))

    for i in 1:nSubject
        tmp_rate[15] = rates[i]
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts[i]]
        K = PBPK(ps)
        u0_eoi = inv(K)*ExponentialAction.expv(durs[i], K, tmp_rate) - inv(K)*tmp_rate

        ut[:,starts[i]:ends[i]] = hcat(ExponentialAction.expv_sequence(times[i], K, u0_eoi)...)

        tmp_sol = ut[15,starts[i]:ends[i]] ./ (VVBs[i]*BP/1000.0)
        append!(predicted, tmp_sol)
    end

    # update that one record that took place before eoi
    tmp_rate_1 = zeros(ncmt)
    tmp_rate_1[15] = rates[1]
    ut_tmp = inv(K_tmp)*ExponentialAction.expv(times[1][1], K_tmp, tmp_rate_1) - inv(K_tmp)*tmp_rate_1
    sol_tmp = ut_tmp[15,1] ./ (VVBs[1]*BP/1000.0)
    predicted2 = [sol_tmp; predicted[2:end]]

    # likelihood
    for i = 1:length(predicted2)
        data[i] ~ LogNormal(log.(max(predicted2[i], 1e-12)), σ)
    end
end

mod = fitPBPK(dat_obs.DV, nSubject, rates, times, durs, wts, VVBs, BP, starts, ends, ncmt)

# run 
nsampl = 250
nchains = 4
adapt_delta = .8
@time mcmcchains = sample(mod, NUTS(nsampl,adapt_delta), MCMCThreads(), nsampl, nchains)  # parallel
@time mcmcchains_prior = sample(mod, Prior(), MCMCThreads(), nsampl, nchains)  # parallel

## save mcmcchains
write(joinpath(modPath, string(modName, "chains.jls")), mcmcchains)
write(joinpath(modPath, string(modName, "chains_prior.jls")), mcmcchains_prior)

##load saved chains
#mcmcchains = read(joinpath(modPath, string(modName, "chains.jls")), Chains)


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
mod_pred = fitPBPK(dat_missing, nSubject, rates, times, durs, wts, VVBs, BP, starts, ends, ncmt)
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

    tmp_rate = zeros(ncmt)
    tmp_rate[15] = rates[1]
    ps_tmp = [CLintᵢ[1], KbBR, KbMU, KbAD, KbBO, KbRB, wts[1]]
    K_tmp = PBPK(ps_tmp)
    u0_eoi_tmp = inv(K_tmp)*ExponentialAction.expv(durs[1], K_tmp, tmp_rate) - inv(K_tmp)*tmp_rate
    ut = zeros(Base.promote_eltype(times[1], K_tmp, u0_eoi_tmp), length(u0_eoi_tmp), length(dat_obs.DV))

    predicted = []

    for i in 1:nSubject
        tmp_rate[15] = rates[i]
        ps = [CLintᵢ[i], KbBR, KbMU, KbAD, KbBO, KbRB, wts[i]]
        K = PBPK(ps)
        u0_eoi = inv(K)*ExponentialAction.expv(durs[i], K, tmp_rate) - inv(K)*tmp_rate
        ut[:,starts[i]:ends[i]] = hcat(ExponentialAction.expv_sequence(times[i], K, u0_eoi)...)

        tmp_sol = ut[15,starts[i]:ends[i]] ./ (VVBs[i]*BP/1000.0)
        append!(predicted, tmp_sol)

        # update that one record that took place before eoi
        tmp_rate_1 = zeros(ncmt)
        tmp_rate_1[15] = rates[1]
        ut_tmp = inv(K_tmp)*ExponentialAction.expv(times[1][1], K_tmp, tmp_rate_1) - inv(K_tmp)*tmp_rate_1
        sol_tmp = ut_tmp[15,1] ./ (VVBs[1]*BP/1000.0)
        predicted2 = [sol_tmp; predicted[2:end]]
    end

    array_pred[j, :] = rand.(LogNormal.(log.(predicted2), df_params[j,:σ]))
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
CSV.write(joinpath(modPath, "df_pred.csv"), df_pred_new2)

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
