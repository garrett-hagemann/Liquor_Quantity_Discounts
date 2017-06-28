module LiquorQD

using Optim

mutable struc Liquor # should be thought of as a product. Products have characteristics
  # maybe should add name as well
  id::Float64
  price::Float64 # price of liquor. Can be changed
  obs_price::Float64 # field so we can reset price to observed
  imported::Bool # imported flag. Should not change
  proof::Float64 # Proof. Should not change
  group_util::Float64 # Utiliy from being in group. i.e. value of d_gin, etc.
  prod_util::Float64 # Product specific utility. i.e. value of d_j
end

mutable struc Market # A market is a collection of liquors
  products::Array{Liquor,1}
  month::Int64
  year::Int64
end

mutable struc DemandCoefs # type to hold estimated demand coefficients
  price::Float64
  imported::Float64
  proof::Float64
end

mutable struc PriceSched # type to hold a price schedule
  rho::Array{Float64,1} # price options
  q_cut::Array{Float64,1} # Quantity cutoffs
  t_cut::Array{Float64,1} # Type cutoffs (lambdas)
  A::Array{Float64,1} # Fixed fees
  N::Int64 # Parts in tariff. An N-part tariff has N-1 prices, FF, etc.
end

function ks_dist(x::Float64,a::Float64,b::Float64)
  if a <= 0
    error("Must have a > 0")
  end
  if b <= 0
    error("Must have b > 0")
  end
  if x < 0 | x > 1
    error("Distribution is defined only on x=[0,1]")
  end
  return 1 - (1-x^a)^b
end

function ks_dens(x::Float64,a::Float64,b::Float64)
  if a <= 0
    error("Must have a > 0")
  end
  if b <= 0
    error("Must have b > 0")
  end
  if x < 0 | x > 1
    error("Density is defined only on x=[0,1]")
  end
  return a*b*x^(a-1)*(1-x^a)^(b-1)
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
  for j in other_products
    res = res + exp(coefs.price*j.price + coefs.imported*j.imported +
      coefs.proof*j.proof + j.group_util + j.prod_util)
  end
  denom = 1 + num + res
  return num/denom
end

function share(price::Float64,product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
  #= coefs argument should be a Lx1 vector where each row contains a pointer to
  a set of coefficients (a DemandCOefs object). The wight vector must be in the
  corresponding order. The length can be as long as needed (i.e. however many
  levels are needed) =#

  product.price = price # changing price of product of interest
  res = 0.0
  for b,w in zip(coefs,weights)
    # b[1] should refer to a DemandCoefs object and b[2] the relevant weight
    res = res + logit_prob(product,b,mkt)*w
  end
  return res
end

function p_star(mc::Float64,product::Liquor,coefs::Array{DemandCoefs,1},weights::Array{Float64,1},mkt::Market)
    # mc = marginal cost. should be in terms of bottles because that's how demand is specified
    function g(p::Float64)
      return (p - mc)*share(p,product,coefs,weights,mkt)
    end
    p_star_optim = Optim.optimize(g,0.0,400.0,method=Brent(), show_trace = false, extended_trace = false, rel_tol = 1e-6, abs_tol = 1e-6)
    return Optim.minimizer(p_star_optim)[1]  # gets to scalar instead of vector.
end



# module end
end
