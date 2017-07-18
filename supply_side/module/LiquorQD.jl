module LiquorQD

using Optim

mutable struct PriceSched # type to hold a price schedule
  rhos::Array{Float64,1} # price options
  t_cuts::Array{Float64,1} # Type cutoffs (lambdas)
  N::Int64 # number of options. A N part tariff as N-1 prices
end

mutable struct Liquor # should be thought of as a product. Products have characteristics
  # maybe should add name as well
  id::Int64 # product number
  price::Float64 # price of liquor. Can be changed
  obs_price::Float64 # field so we can reset price to observed
  imported::Bool # imported flag. Should not change
  proof::Float64 # Proof. Should not change
  group_util::Float64 # Utiliy from being in group. i.e. value of d_gin, etc.
  prod_util::Float64 # Product specific utility. i.e. value of d_j
  ps::Nullable{PriceSched} # holds price schedule for this product. Can be null if no sched matched
end

mutable struct Market # A market is a collection of liquors
  products::Array{Liquor,1}
  month::Int64
  year::Int64
end

mutable struct DemandCoefs # type to hold estimated demand coefficients
  price::Float64
  imported::Float64
  proof::Float64
end

mutable struct WholesaleParams # type to hold wholesaler parameters
  c::Float64 # marginal cost for wholesaler
  a::Float64 # a parameter for the Kumaraswamy distribution faced by wholesaler
  b::Float64 # b parameter for the Kumaraswamy distribution faced by wholesaler
end

function ks_dist(x::Float64,a::Float64,b::Float64)
  if a <= 0.0
    error("Must have a > 0")
  end
  if b <= 0.0
    error("Must have b > 0")
  end
  if (x < 0.0)
    return 0.0
  end
  if (x > 1.0)
    return 1.0
  end
  return 1.0 - (1.0-x^a)^b
end

function ks_dens(x::Float64,a::Float64,b::Float64)
  if a <= 0.0
    error("Must have a > 0")
  end
  if b <= 0.0
    error("Must have b > 0")
  end
  if (x < 0.0) | (x > 1.0)
    #error("Density is defined only on x=[0,1]. Tried to evaluate at $x")
    return 0.0
  end
  return a*b*x^(a-1.0)*(1.0-x^a)^(b-1.0)
end

function sparse_int(f::Function, a::Float64, b::Float64)
	#= Implements sparse grid quadrature from sparse-grids.de
	This implements the 1 dimensional rule that is exact for
	polynomials up to order 11 (using the accuracy level 6 rule).
  This version provides comparable results to the built in quadgk, but is
  slightly faster due to fewer function calls and no calculation of
  approximation error. Deviations in accuracy and speed are minor at best.

	f = one dimensional function
	a = lower bound of integration
	b = upperbound of integration	=#
	weights = [.052328105232810528, .13424401342440137, .20069872006987202, .22545832254583226, .20069872006987202, .13424401342440137, .052328105232810528]
	nodes = [.01975439999999995,.11270170000000002,.28287810000000002,.5,.71712189999999998, .88729829999999998,.98024560000000005]
	return dot([f((b-a)*u + a) for u in nodes], weights)*(b-a)::Float64
end

function logit_prob(product::Liquor, coefs::DemandCoefs, mkt::Market)
  #= The estimated share function has some discrete variation in addition to
    possible random coefficients. In particular, the price, proof, and imported
    coefficients vary according to 4 income categories. Thus, for a given
    set of product characteristics, there are 4 evaluations of the logit share
    function that must be evaluated then weighted. This function does the base
    logit share function evaluation. the actual share function will then be a
    weighted evaluation of these. The random coefficient version can then be
    an integral of that weighted sum over the distribution of coefficients
    for each type.
  =#
  other_prods = filter(e->e!=product,mkt.products)
  num = exp(coefs.price*product.price + coefs.imported*product.imported +
    coefs.proof*product.proof + product.group_util + product.prod_util)

  res = 0.0
  for j in other_prods
    res = res + exp(coefs.price*j.price + coefs.imported*j.imported +
      coefs.proof*j.proof + j.group_util + j.prod_util)
  end
  denom = 1 + num + res
  return num/denom
end

function share(price::Float64,product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
  #= coefs argument should be a Lx1 vector where each row contains a pointer to
  a set of coefficients (a DemandCoefs object). The wight vector must be in the
  corresponding order. The length can be as long as needed (i.e. however many
  levels are needed) =#

  product.price = price # changing price of product of interest
  res = 0.0
  for (b,w) in zip(coefs,weights)
    # b[1] should refer to a DemandCoefs object and b[2] the relevant weight
    res = res + logit_prob(product,b,mkt)*w
  end
  return res
  product.price = product.obs_price # undoing change to price for subsequent calls
end

function d_share(price::Float64,product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
  product.price = price
  res = 0.0
  for (b,w) in zip(coefs,weights)
    s = logit_prob(product,b,mkt)
    res = res + b.price*s*(1-s)*w
  end
  return res
  product.price = product.obs_price
end

function p_star(mc::Float64,product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
    # mc = marginal cost. should be in terms of bottles because that's how demand is specified
    #=
    # Optimizer approach
    function g(p::Float64)
      return -(p - mc)*share(p,product,coefs,weights,mkt) # want to return neg of profit because we are going to minimize it
    end
    p_star_optim = Optim.optimize(g,0.0,10.0*mc,method=Brent(), show_trace = false, extended_trace = false, rel_tol = 1e-6, abs_tol = 1e-6)
    return Optim.minimizer(p_star_optim)[1]  # gets to scalar instead of vector.
    =#

    # Fixed point approach
    g(p::Float64) = -share(p,product,coefs,weights,mkt)/d_share(p,product,coefs,weights,mkt) + mc
    diff = 1.0
    tol = 1e-12 # want tight tolerance so we don't popigate numerical errors
    old_p = mc
    while diff > tol
      new_p = g(old_p)
      diff = abs(new_p - old_p)
      old_p = new_p
    end
    return old_p

end

function obs_lambdas(rhos::Array{Float64,1},ff::Array{Float64,1},product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market,M::Float64)
  # function calculates type cutoffs from observed FF and prices
  obs_t = [0.0] # corresponds to lambda1. Must be zero as A0 = 0 & A1 = 0
  ps(mc::Float64) = p_star(mc,product,coefs,weights,mkt)
  s(p::Float64) = share(p, product,coefs,weights,mkt)
  for k = 2:length(rhos)
    res = (ff[k] - ff[k-1])/(M*((ps(rhos[k]) - rhos[k])*s(ps(rhos[k])) - (ps(rhos[k-1]) - rhos[k-1])*s(ps(rhos[k-1]))))
    push!(obs_t,res)
  end
  return obs_t
end

function wholesaler_profit(ps::PriceSched, params::WholesaleParams, product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market,ps_pre::Dict{Int64,Float64}=Dict{Int64,Float64}())
  #= ps_pre holds pre-calculated values of p_star. This optional argument prevents recalculation as the p_star values
  don't depend on any wholesaler parameters =#
  
  M = 10000.0 # normalizing constant. Shouldn't affect price schedule, just scale of profits
  lambda_vec = [ps.t_cuts; 1.0] # need to put in upper boundary values
  rho_vec = ps.rhos # no boundary here, will deal with in iterator
  N = ps.N
  # defining some functions for convience
  s(p) = share(p,product,coefs,weights,mkt)
  est_pdf(x) = ks_dens(x,params.a,params.b)
  est_cdf(x) = ks_dist(x,params.a,params.b)
  profit = 0.0
  for i = 1:N-1
    k = i # dealing with indexing
    f(l) =  l*est_pdf(l)
    int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
    if (k == 1)
      if isempty(ps_pre)
        ps1 = p_star(rho_vec[k],product,coefs,weights,mkt)
      else
        ps1 = ps_pre[k]
      end
      inc = ((rho_vec[k] - params.c)*M*s(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*M*((ps1 - rho_vec[k])*s(ps1) - 0.0)
      profit = profit + inc
    else
      # Pre-calculating some stuff to avoid repeated calls to p_star
      if isempty(ps_pre)
        ps1 = p_star(rho_vec[k],product,coefs,weights,mkt)
        ps2 = p_star(rho_vec[k-1],product,coefs,weights,mkt)
      else
        ps1 = ps_pre[k]
        ps2 = ps_pre[k-1]
      end
      inc = ((rho_vec[k] - params.c)*M*s(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*M*((ps1 - rho_vec[k])*s(ps1) - (ps2 - rho_vec[k-1])*s(ps2))
      profit = profit + inc
    end
  end
  return profit
  #=
  #### Linear version. Only for checking with Wilson Problem
  M = 1.0 # normalizing constant. Shouldn't affect price schedule, just scale of profits
  lambda_vec = [ps.t_cuts; 1.0] # need to put in upper boundary values
  rho_vec = ps.rhos # no boundary here, will deal with in iterator
  N = ps.N
  # defining some functions for convience
  s(p) = max((1-p),0.0)
  est_pdf(x) = ks_dens(x,params.a,params.b)
  est_cdf(x) = ks_dist(x,params.a,params.b)
  profit = 0.0
  for i = 1:N-1
    k = i # dealing with indexing
    f(l) =  l*est_pdf(l)
    int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
    if (k == 1)
      ps1 = (1 + (rho_vec[k]))/2
      inc = ((rho_vec[k] - params.c)*M*s(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*M*((ps1 - rho_vec[k])*s(ps1) - 0.0)
      profit = profit + inc
    else
      # Pre-calculating some stuff to avoid repeated calls to p_star
      ps1 = (1 + (rho_vec[k]))/2
      ps2 = (1 + (rho_vec[k-1]))/2
      inc = ((rho_vec[k] - params.c)*M*s(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*M*((ps1 - rho_vec[k])*s(ps1) - (ps2 - rho_vec[k-1])*s(ps2))
      profit = profit + inc
    end
  end
  return profit
  =#
end

function optimal_price_sched(params::WholesaleParams, N::Int64, product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
  #= params contains the parameters for the wholesaler. This means the marginal
  cost as well as the parameters of the Kumaraswamy distribution. =#
  function g(x::Array{Float64,1})
    ps = PriceSched(x[1:N-1],[0.0; x[N:end]],N) # intializing new price schedule object. Need to impose constrained
    if !all((0.0 .<= ps.t_cuts .<= 1.0)) | !all((0.0 .<= ps.rhos)) # error checking to limit derivative-free search
  		return Inf
  	else
      return -1.0*wholesaler_profit(ps,params,product,coefs,weights,mkt) # return wholesaler profit for price schedule. negative because use minimizer
    end
  end
  rho_guess = sort(rand(N-1),rev=true) # prices in 0,1 declining. Scaling up seems to have no effect on solutions (which is good)
  t_guess = sort(rand(N-2)) # generates types between 0 and 1 that increase. One less b/c of constrained
  hs_x0 = [rho_guess; t_guess]
  hs_res = Optim.optimize(g,hs_x0,method=SimulatedAnnealing(),iterations=1) # taking 25000 SA steps to get good guess. stabalizes the solution
  x0 = Optim.minimizer(hs_res)
  ps_optim_res = Optim.optimize(g,x0,method=NelderMead())
  optim_rho = Optim.minimizer(ps_optim_res)[1:N-1]
  optim_t = Optim.minimizer(ps_optim_res)[N:end]
  return PriceSched(optim_rho,[0.0; optim_t],N) # constrained
end

function dev_gen(ps::PriceSched,δ::Float64)
  #= returns deviations of a price schedule. The deviations constitute all price
  schedules where one or more elements is multiplied by 1+δ or 1-δ. Function
  should output array of price schedule objects =#
  N = ps.N
  dev_steps = [1+δ,1.0,1-δ] # deviation amounts
  dev_poss = fill(dev_steps,2*(N-1)) # deviation possibilities for each ps element
  ps_vec = [ps.rhos ; ps.t_cuts] # price schedule in vector form so we can manipulate it
  dev_cart_prod = Iterators.product(dev_poss...) # cartesian product of elems in dev_poss. Each is a possible deviation vector. Iterators is in base
  dev_cart_prod_filter = Iterators.filter(x->(x!=tuple(ones(2*(N-1))...)),dev_cart_prod) # removing deviation vector of all ones (i.e. unperturbed price schedule)
  dev_ps = map(x->ps_vec.*x,dev_cart_prod_filter) # generating perturbed price schedules
  output = PriceSched[] # initializing empty arry of price schedules.
  Threads.@threads for s in dev_ps  # potentially very many deviations. Scales at the rate of (2*(N-1))^3 - 1. Threads may be faster
    push!(output,PriceSched(s[1:N-1],s[N:end],N))
  end
  return output::Array{PriceSched,1}
end

function moment_obj_func(ps::PriceSched, devs::Array{PriceSched,1},params::WholesaleParams,product::Liquor,coefs::Array{DemandCoefs,1},weight::Array{Float64,1},mkt::Market,ps_pre_array::Array{Dict{Int64,Float64}}=Dict{Int64,Float64}[])
  # calculates the value of the objective function for a given set of deviations.
  wp(x::PriceSched,ps_pre::Dict{Int64,Float64}=Dict{Int64,Float64}()) = wholesaler_profit(x,params,product,coefs,weight,mkt,ps_pre)
  obs_profit = wp(ps)
  if obs_profit < 0.0
    res = Inf
  else
    #dev_profit = max.(map(wp,devs),0.0) # likewise with all deviations. Use max.() because array
    dev_profit = map(wp,devs,ps_pre_array) # likewise with all deviations. Use max.() because array
    nu = min.(obs_profit - dev_profit,0.0).^2 # if dev profit is greater, then you get a positive value for that deviation.
    res = sum(nu)
  end
  return res
end

function optimize_moment(ps::PriceSched, devs::Array{PriceSched,1},product::Liquor,coefs::Array{DemandCoefs,1},weight::Array{Float64,1},mkt::Market,iters::Int64,ps_pre_array::Array{Dict{Int64,Float64}}=Dict{Int64,Float64}[];x0=[0.0,1.0,1.0])

  function Q(x::Array{Float64,1})
    theta = WholesaleParams(x...) # search all 3 params
    if (theta.b <= 0.0 ) | (theta.b >= 10.0) | (theta.c < 0.0) | (theta.a < 1.0) # constraining SA search
      return Inf
    else
      return moment_obj_func(ps,devs,theta,product,coefs,weight,mkt,ps_pre_array)
    end
  end

  moment_res = Optim.optimize(Q,x0,method=SimulatedAnnealing(), store_trace=true, show_trace = true, extended_trace=true, iterations=iters)
  moment_res_min = Optim.minimizer(moment_res)
  x_trace = Optim.x_trace(moment_res)
  f_trace = Optim.f_trace(moment_res)
  return (WholesaleParams(moment_res_min...),x_trace,f_trace)
  #=
  function Q(x::Float64)
    theta = WholesaleParams(x,1.0,1.0) # search just c
    if (theta.b <= 0.0 ) | (theta.b >= 30.0) | (theta.c < 0.0) | (theta.a < 1.0) # constraining SA search
      return Inf
    else
      return moment_obj_func(ps,devs,theta,product,coefs,weight,mkt)
    end
  end
  x0 = 0.0
  moment_res = Optim.optimize(Q,0.0,10.0, store_trace=true, show_trace = true, iterations=iters)
  moment_res_min = Optim.minimizer(moment_res)
  return WholesaleParams(moment_res_min,1.0,1.0)
  =#
end




# Export statements

export Liquor, Market, DemandCoefs, WholesaleParams, PriceSched
export ks_dist, ks_dens, sparse_int, share, d_share,
  p_star, wholesaler_profit, optimal_price_sched, dev_gen, moment_obj_func,
  optimize_moment, obs_lambdas

end
