@everywhere using DataArrays, DataFrames, ForwardDiff, Distributions, Optim, Iterators 
#= goofy conversion magic =#
import Base.convert
function convert(t::Type{Float64},x::Array{Float64,1})
	return convert(Float64,x[1])
end

import Base.isless
function isless(x::Array{Float64,1},y::Float64)
	return isless(x[1],y)
end
function isless(x::Array{Float64,1},y::Array{Float64,1})
	return isless(x[1],y[1])
end

srand(69510606) #seeding random number gen
#=
df = readtable("../../demand_estimation/berry_logit/berry_logit.csv")
char_list = [:price, :d_gin, :d_vod, :d_rum, :d_sch, :d_brb, :d_whs, :holiday]
coefs = readtable("../../demand_estimation/berry_logit/blogit_coefs.csv")

#= each row in df represents one bundle of p,q. Each row also contains all
the needed information to determine the cutoffs. So, we'll be looping over the relevant
price schedule segments and cutoffs and doing it for each row. =#

#= each column in coefs is the appropriate coefficient =#
coef_vec = convert(Array,coefs) # row vector
alpha = coef_vec[1,1]
=#
function weibull_pdf(x, gamma, mu,a)
	return (gamma/a)*((x - mu)/a)^(gamma - 1)*exp(-((x - mu)/a)^gamma)
end

function unif_pdf(x)
	if 0 <= x <= 1
		return 1
	else
		return 0
	end
end

function burr_cdf(x,gamma)
	return 1 - (1-x)^(1/gamma)
end

@everywhere function sparse_int(f::Function, a::Float64, b::Float64)
	#= Implements sparse grid quadrature from sparsegrids.de
	This implements the 1 dimensional rule that is exact for
	polynomials up to order 25.

	f = one dimensional function
	a = lower bound of integration
	b = upperbound of integration

	=#
	weights = [.052328105232810528, .13424401342440137, .20069872006987202, .22545832254583226, .20069872006987202, .13424401342440137, .052328105232810528]

	nodes = [.01975439999999995,.11270170000000002,.28287810000000002,.5,.71712189999999998, .88729829999999998,.98024560000000005]

	#f_evals = [f((b-a)*u + a) for u in nodes]
	return dot([f((b-a)*u + a) for u in nodes], weights)*(b-a)::Float64
end

@everywhere function powerset(a)
	res = Any[]
	for i = 1:length(a)
		for j in combinations(a,i)
			push!(res,j)
		end
	end
	return res
end
		
@everywhere function ss_est(mkt_prod_tup::Tuple)
	market = mkt_prod_tup[1]
	product = mkt_prod_tup[2]
	println("Working with Market: $market, Product: $product")
	# defining selection boolians so that we can select appropriate products/markets
	df = readtable("../../demand_estimation/berry_logit/berry_logit.csv")
	char_list = [:price, :d_gin, :d_vod, :d_rum, :d_sch, :d_brb, :d_whs, :holiday]
	coefs = readtable("../../demand_estimation/berry_logit/blogit_coefs.csv")

	#= each row in df represents one bundle of p,q. Each row also contains all
	the needed information to determine the cutoffs. So, we'll be looping over the relevant
	price schedule segments and cutoffs and doing it for each row. =#

	#= each column in coefs is the appropriate coefficient =#
	coef_vec = convert(Array,coefs) # row vector
	alpha = coef_vec[1,1]
	prod_bool = ((df[:product] .== product) & (df[:mkt] .== market))
	other_bool = ((df[:product] .!= product) & (df[:mkt] .== market))

	prod_chars = convert(Array,df[prod_bool,char_list])
	other_chars = convert(Array,df[other_bool,char_list])

	# flag that determines if product has matched price data
	matched_upc = (df[prod_bool, :_merge_purchases][1] == 3)

	if matched_upc == true # Only need to work with matched products. No price sched for non-matched.
		println("Product has matched price data.")
		#Defining some constants for use in the share function
		xb_prod = (prod_chars[:,2:end]*coef_vec[:,2:end]')[1] + df[prod_bool,:xi][1] # scalar
		
		if size(other_chars)[1] == 0
			xb_others = 0
		else
			xb_others = other_chars*coef_vec' + convert(Vector,df[other_bool,:xi]) # J-1 x 1
		end

		#Grabbing market size. Needs to be the same as that used to estimate shares
		M = df[prod_bool,:M][1]
		#M = 3000
		tol = 1e-16*M # tolerance for price schedule optimization solution
		inner_tol = 1e-2 # tolerance for price schedule optimization solution
		prod_price = df[prod_bool, :price][1]

		# lists of observed prices and quantities. Need these to calculate error using 
		# estimated price schedule 
		rho_list = [:actual_p0, :actual_p1, :actual_p2, :actual_p3, :actual_p4, :actual_p5, :actual_p6,
			:actual_p7, :actual_p8, :actual_p9, :actual_p10]

		cutoff_q_list = [:disc_q0, :disc_q1, :disc_q2, :disc_q3, :disc_q4, :disc_q5, :disc_q6,
			:disc_q7, :disc_q8, :disc_q9, :disc_q10]

		obs_rhos = dropna(convert(DataArray,df[prod_bool, rho_list])'[:,1]) #some goofy conversion to make the dropna work as expected
		obs_cutoff_q = dropna(convert(DataArray,df[prod_bool, cutoff_q_list])'[:,1]) #same as above
		obs_ff = [0.0] # first fixed fee is by definition, and corresponds to A1
		
		max_rho = copy(maximum(obs_rhos))
		max_mc = 15.0
		# calculating series of A (intercepts of piecewise linear price schedules
		for k = 2:length(obs_cutoff_q)
			diff = (obs_cutoff_q[k] - 1)*(obs_rhos[k-1] - obs_rhos[k])
			res = diff + obs_ff[k-1]
			push!(obs_ff,res)
		end


		#Defining share function which is ONLY a function of price
		function share(p)
			#= Function calculates shares of a given  product in a given market 
			based on price. Characteristics of the product and all other products
			are assumed fixed. 

			p = price =#
			
			num = exp(alpha*p + xb_prod)
			denom = 1 + num + sum(exp(xb_others))
			s = num/denom	

			return s[1] # dealing with weird type promotion stuff
		end
		#Derivatives of the share fuction
		d_share(p) = ForwardDiff.derivative(share, p[1]) #p should be scalar. Optimization routine returns a 1 element array at some point
		dd_share(p) = ForwardDiff.derivative(d_share,p[1])
		ddd_share(p) = ForwardDiff.derivative(dd_share,p[1])

		### Optimizer pstar
		function p_star(rho)
			function g(p::Float64) # retailer profit. Note that optimal p is independent of M and type
				return -((p - rho)*share(p))[1] # Ignores fixed cost as it just shifts problem. Need the indexing for some reason in optimize
			end
			function gj!(p,gvec)
				gvec[1] =  -((p - rho)*d_share(p) + share(p))
			end
			function gh!(p, gjvec)
				gjvec[1] = -((p - rho)*dd_share(p) + 2.0*d_share(p))
			end
			poptim_res = Optim.optimize(g,0.0,100.0, method=Brent(), show_trace = false, extended_trace = false, rel_tol = 1e-6, abs_tol = 1e-6)
			#poptim_res = optimize(uLg,focs!,ux0,BFGS(),OptimizationOptions(show_every = true, extended_trace = true, iterations = 1500, g_tol = 1e-6))
			return Optim.minimizer(poptim_res)[1]
		end
		function d_pstar_d_rho(rho) 
			u = p_star(rho)
			res = d_share(u) ./ (dd_share(u)*(u - rho) + 2*d_share(u))
			return res[1]
		end
		function d2_pstar_d2_rho(rho)
			u = p_star(rho)
			num1 = (dd_share(u)*(u - rho) + 2*d_share(u))*dd_share(u)*d_pstar_d_rho(rho)
			num2 = d_share(u)*(dd_share(u)*(d_pstar_d_rho(rho) - 1) + ((u - rho)*ddd_share(u) + 2*dd_share(u))*d_pstar_d_rho(rho))
			den = dd_share(u)*(u - rho) + 2*d_share(u)
			res = (num1 - num2)/(den^2)
			return res
		end
		
		# Given observed As and observed rhos, can calculate implied lambdas
		obs_lambdas = [0.0] # corresponds to lambda1. Must be zero as A0 = 0 & A1 = 0
		for k = 2:length(obs_rhos)
			res = (obs_ff[k] - obs_ff[k-1])/(M*((p_star(obs_rhos[k]) - obs_rhos[k])*share(p_star(obs_rhos[k])) - (p_star(obs_rhos[k-1]) - obs_rhos[k-1])*share(p_star(obs_rhos[k-1]))))
			push!(obs_lambdas,res)
		end
		obs_sched = [obs_rhos;obs_lambdas]
		println(obs_sched)
		#M = 1

		# Defining parameters for linear approx to demand to give hot start to non-linear version.
		# apporximating demand around observed product price
		LA = share(prod_price) - d_share(prod_price)*prod_price
		LB = d_share(prod_price)
		function Lshare(p)
			return LA + LB*p
		end
		Ld_share(p) = LB
		Ldd_share(p) = 0.0
		Lddd_share(p) = 0.0

		function Lp_star(rho)
			return -LA/(2.0*LB) + (rho)/2.0
		end
		function Ld_pstar_d_rho(rho)
			return 0.5
		end
		function Ld2_pstar_d2_rho(rho)
			return 0.0
		end

		function Lprice_sched_calc(params::Array{Float64,1},N::Int)
			constrained = 1 # solving problem with A1 constrained to 0
			m = 3000

			# params are the params of the wholesaler's problem we're trying to
			# estimate. 

			# Function returns the coefficients defining the price schedule as [rho, lambda]
			# where rho is n long and lambda is n-1 long for an n option schedule

			c = params[1] # marginal cost for wholesaler
			#max_mc = params[2] # max mC for retailer. Scales type parameter
			#m = exp(params[5])
			#distribution parameters
			a = exp(params[2]) #
			b = exp(params[3]) # 
			lambda_lb = 0.0
			lambda_ub = 1.0
				
			est_cdf(x) = cdf(Beta(a,b),x)
			est_pdf(x) = pdf(Beta(a,b),x)
			d_est_pdf(x) = est_pdf(x)*((a + b - 2)*x - (a-1))/((x-1)*x)
		
			if (constrained == 1) # contrained => A1 = 0 <=> lambda_1 = 0
				function urho(k::Integer)
					return -LA/LB + (LA + LB*c)/(2*LB)*(2*k - 1)/(N-1)
				end
				function ulambda(k::Integer)
					return (k-1)/(N-1)
				end
				rho_0 = urho(0)
				function Lw_profit(theta) #theta should contain rho_1, ..., rho_n-1 ; lambda_2, ... , lambda_n-1
					lambda_vec = [lambda_lb; 0.0; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]]
					profit = 0.0
					for i = 1:N-1
						k = i+1 # dealing with indexing
						f(l) =  l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						inc = ((rho_vec[k] - c)*m*Lshare(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*Lshare(ps1)*m - (ps2 - rho_vec[k-1])*Lshare(ps2)*m)
						profit = profit + inc
					end
					return -profit
				end
				function Lwfocs!(theta, wfocs_vec)
					lambda_vec = [lambda_lb; 0.0 ; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					# Calculating FOCs
					for i in 1:N-1 
						k = i+1 # dealing with the increment in indexing and making code look more like notation
						#calculating integral for rho FOC
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						#rho FOC
							term1 = ((rho_vec[k] - c)*m*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + m*Lshare(ps1))*int
							term2 = ((ps1 - rho_vec[k])*m*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + m*Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
							term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
							res = term1 + term2*term3
							wfocs_vec[i] = -res
						# lambda FOC
						if i == 1
							# do nothing
						else
							term1 = ((rho_vec[k] - c)*m*Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((rho_vec[k-1] - c)*m*Lshare(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
							term3 = ((ps1 - rho_vec[k])*m*Lshare(ps1) - (ps2 - rho_vec[k-1])*m*Lshare(ps2))
							term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2 + term3*term4
							wfocs_vec[i+N-1-1] = -res
						end
						
					end
				end
				function Lwhess!(theta, whess_mat)
					lambda_vec = [lambda_lb; 0.0; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					
					#rho rho part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*(Ld_share(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + Ld_pstar_d_rho(rho_vec[k])^2*Ldd_share(ps1)) + 2*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*(Ld_share(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + Ld_pstar_d_rho(rho_vec[k])^2*Ldd_share(ps1)) + Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k])*(Ld_pstar_d_rho(rho_vec[k]) - 1) + Lshare(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + (Ld_pstar_d_rho(rho_vec[k]) - 1)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]))
								term3 = (1-est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
								res = term1*int + term2*term3
								whess_mat[i,j] = -res
							else
								whess_mat[i,j] = 0.0
							end
						end
					end
					# rho lambda part of hess (upper right)
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 2:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k])) + lambda_vec[k]*(-est_pdf(lambda_vec[k]))
								res = term1 + term2*term3
								whess_mat[i,j+N-1-1] = -res
							elseif j == i+1
								term1 = m*((rho_vec[k] - c)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1))*(lambda_vec[k+1]*est_pdf(lambda_vec[k+1]))
								term2 = m*((ps1 - rho_vec[k])*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k+1])) + lambda_vec[k+1]*(-est_pdf(lambda_vec[k+1]))
								res = term1 + term2*(-1)*term3
								whess_mat[i,j+N-1-1] = -res
									
							else
								whess_mat[i,j+N-1-1] = 0.0
							end
						end
					end
					#lambda lambda part of hess
					for i = 2:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 2:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*Lshare(ps1))
								term2 = (-(lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k])))
								term3 = m*((rho_vec[k-1] - c)*Lshare(ps2))
								term4 = (lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k]))
								term5 = m*((ps1 - rho_vec[k])*Lshare(ps1) - (ps2 - rho_vec[k-1])*Lshare(ps2))
								term6 = -est_pdf(lambda_vec[k]) + lambda_vec[k]*(-d_est_pdf(lambda_vec[k])) + (-est_pdf(lambda_vec[k]))
								res = term1*term2 + term3*term4 + term5*term6
								whess_mat[i+N-1-1,j+N-1-1] = -res
							else
								whess_mat[i+N-1-1,j+N-1-1] = 0.0
							end
						end
					end
					# lambda rho part (lower left). should be same as rho_lambda part.
					for i = 2:N-1
						for j = 1:N-1
							whess_mat[i+N-1-1,j] = whess_mat[j,i+N-1-1]
						end
					end
				end
			else
				function urho(k::Integer)
					return -LA/LB + (1/2*N-1)*(2/LB)*(LA + LB*c)*k
				end
				function ulambda(k::Integer)
					return (2*k - 1)/(2*N-1)
				end
				rho_0 = urho(0)
				function Lw_profit(theta) # params should contain rho_1 ... rho_N-1, lambda_1, ..., lambda n-1
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]]
					profit = 0.0
					for i = 1:N-1
						k = i+1 # dealing with indexing
						f(l) =  l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						inc = ((rho_vec[k] - c)*m*Lshare(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*Lshare(ps1)*m - (ps2 - rho_vec[k-1])*Lshare(ps2)*m)
						profit = profit + inc
					end
					return -profit
				end
				function Lwfocs!(theta, wfocs_vec)
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					# Calculating FOCs
					for i in 1:N-1 
						k = i+1 # dealing with the increment in indexing and making code look more like notation
						#calculating integral for rho FOC
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						#rho FOC
							term1 = ((rho_vec[k] - c)*m*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + m*Lshare(ps1))*int
							term2 = ((ps1 - rho_vec[k])*m*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + m*Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
							term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
							res = term1 + term2*term3
							wfocs_vec[i] = -res
						# lambda FOC
						if i == 1
							term1 = ((rho_vec[k] - c)*m*Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((ps1 - rho_vec[k])*m*Lshare(ps1) - (ps2 - rho_vec[k-1])*m*Lshare(ps2))
							term3 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2*term3
							wfocs_vec[i+N-1] = -res
						else
							term1 = ((rho_vec[k] - c)*m*Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((rho_vec[k-1] - c)*m*Lshare(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
							term3 = ((ps1 - rho_vec[k])*m*Lshare(ps1) - (ps2 - rho_vec[k-1])*m*Lshare(ps2))
							term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2 + term3*term4
							wfocs_vec[i+N-1] = -res
						end
						
					end
				end
				function Lwhess!(theta,whess_mat)
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					
					#rho rho part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*(Ld_share(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + Ld_pstar_d_rho(rho_vec[k])^2*Ldd_share(ps1)) + 2*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*(Ld_share(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + Ld_pstar_d_rho(rho_vec[k])^2*Ldd_share(ps1)) + Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k])*(Ld_pstar_d_rho(rho_vec[k]) - 1) + Lshare(ps1)*Ld2_pstar_d2_rho(rho_vec[k]) + (Ld_pstar_d_rho(rho_vec[k]) - 1)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]))
								term3 = (1-est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
								res = term1*int + term2*term3
								whess_mat[i,j] = -res
							else
								whess_mat[i,j] = 0.0
							end
						end
					end
					# rho lambda part of hess (upper right)
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k])) + lambda_vec[k]*(-est_pdf(lambda_vec[k]))
								res = term1 + term2*term3
								whess_mat[i,j+N-1] = -res
							elseif j == i+1
								term1 = m*((rho_vec[k] - c)*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1))*(lambda_vec[k+1]*est_pdf(lambda_vec[k+1]))
								term2 = m*((ps1 - rho_vec[k])*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k+1])) + lambda_vec[k+1]*(-est_pdf(lambda_vec[k+1]))
								res = term1 + term2*(-1)*term3
								whess_mat[i,j+N-1] = -res
									
							else
								whess_mat[i,j+N-1] = 0.0
							end
						end
					end
					#lambda lambda part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = Lp_star(rho_vec[k])
						ps2 = Lp_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*Lshare(ps1))
								term2 = (-(lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k])))
								term3 = m*((rho_vec[k-1] - c)*Lshare(ps2))
								term4 = (lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k]))
								term5 = m*((ps1 - rho_vec[k])*Lshare(ps1) - (ps2 - rho_vec[k-1])*Lshare(ps2))
								term6 = -est_pdf(lambda_vec[k]) + lambda_vec[k]*(-d_est_pdf(lambda_vec[k])) + (-est_pdf(lambda_vec[k]))
								res = term1*term2 + term3*term4 + term5*term6
								whess_mat[i+N-1,j+N-1] = -res
							else
								whess_mat[i+N-1,j+N-1] = 0.0
							end
						end
					end
					# lambda rho part (lower left). should be same as rho_lambda part.
					for i = 1:N-1
						for j = 1:N-1
							whess_mat[i+N-1,j] = whess_mat[j,i+N-1]
						end
					end
				end
			end
			##### Optimizer Approach
			rho_start = convert(Array{Float64,1},[urho(k) for k = 1:N-1])
			lambda_start = convert(Array{Float64,1},[ulambda(k) for k = 2:N-1])	
			innerx0 = [rho_start ; lambda_start]
			# checking hessian and gradient
			#=	
			innerx0 = [20.0, 10.0, 5.0, 3.0, 0.3, 0.6, 0.9]
			println(innerx0)
			gtest1 = ones(length(innerx0))
			gtest2 = ones(length(innerx0))
			htest = ones(length(innerx0),length(innerx0))
			eps = zeros(length(innerx0)) 
			step = 1e-9
			eps[7] = step
			upx = innerx0 + eps
			downx = innerx0 - eps
			est_grad = (Lw_profit(upx) - Lw_profit(downx))/(2*step)
			println("Numerical grad: ", est_grad)
			Lwfocs!(innerx0,gtest1)
			println(gtest1)
			
			Lwfocs!(upx,gtest1)
			Lwfocs!(downx,gtest2)
			Lwhess!(innerx0,htest)
			println("Numerical Hessian: ", (gtest1 - gtest2)/(2*step))
			println(htest)
			=#
			res = Optim.optimize(Lw_profit,Lwfocs!,Lwhess!,innerx0,method=NewtonTrustRegion())
			sol_sched=Optim.minimizer(res)
			est_rhos = [rho_0 ; sol_sched[1:N-1]]
			est_lambdas = [lambda_lb ; 0.0; sol_sched[N:end] ; lambda_ub]
			
			# Calculating fixed fees
			est_ff = [0.0]
			for i in 1:N-1
				k = i+1
				A = est_ff[k-1] + (Lp_star(est_rhos[k]) - est_rhos[k])*Lshare(Lp_star(est_rhos[k]))*m*est_lambdas[k] - (Lp_star(est_rhos[k-1]) - est_rhos[k-1])*Lshare(Lp_star(est_rhos[k-1]))*m*est_lambdas[k]
				push!(est_ff,A)
			end
			return (est_rhos,est_ff,est_lambdas)
		end		
		
		function price_sched_calc(params,N)
			constrained = 1 # solving problem with A1 constrained to 0
			m = 3000 # normalizing as it shouldn't matter for schedule calculation. Helps numerical methods apparently
			# params are the params of the wholesaler's problem we're trying to
			# estimate. 

			# Function returns the coefficients defining the price schedule as [rho, lambda]
			# where rho is n long and lambda is n-1 long for an n option schedule

			c = params[1] # marginal cost for wholesaler
			#max_mc = params[2] # max MC for retailer. Scales type parameter
			#M = exp(params[5])
			#distribution parameters
			a = exp(params[2]) #
			b = exp(params[3]) # 
			lambda_lb = 0.0
			lambda_ub = 1.0
				
			est_cdf(x) = cdf(Beta(a,b),x)
			est_pdf(x) = pdf(Beta(a,b),x)
			d_est_pdf(x) = est_pdf(x)*((a + b - 2)*x - (a-1))/((x-1)*x)
			
			rho_0 = 2.0*max_rho
		
			##### Optimizer Approach
			if (constrained == 1)
				function w_profit(theta) #theta should contain rho_1, ..., rho_n-1 ; lambda_2, ... , lambda_n-1
					lambda_vec = [lambda_lb; 0.0 ; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]]
					profit = 0.0
					for i = 1:N-1
						k = i+1 # dealing with indexing
						f(l) =  l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						inc = ((rho_vec[k] - c)*m*share(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*share(ps1)*m - (ps2 - rho_vec[k-1])*share(ps2)*m)
						profit = profit + inc
					end
					return -profit
				end
				function wfocs!(theta, wfocs_vec)
					lambda_vec = [lambda_lb; 0.0; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					# Calculating FOCs
					for i in 1:N-1 
						k = i+1 # dealing with the increment in indexing and making code look more like notation
						#calculating integral for rho FOC
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						#rho FOC
							term1 = ((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*int*m
							term2 = ((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))*m
							term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
							res = term1 + term2*term3
							wfocs_vec[i] = -res
						# lambda FOC
						if i == 1
							# do nothing
						else
							term1 = ((rho_vec[k] - c)*m*share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((rho_vec[k-1] - c)*m*share(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
							term3 = ((ps1 - rho_vec[k])*m*share(ps1) - (ps2 - rho_vec[k-1])*m*share(ps2))
							term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2 + term3*term4
							wfocs_vec[i+N-1-1] = -res
						end
						
					end
				end
				function whess!(theta, whess_mat)
					lambda_vec = [lambda_lb; 0.0 ; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					
					#rho rho part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*(d_share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + d_pstar_d_rho(rho_vec[k])^2*dd_share(ps1)) + 2*d_share(ps1)*d_pstar_d_rho(rho_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*(d_share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + d_pstar_d_rho(rho_vec[k])^2*dd_share(ps1)) + d_share(ps1)*d_pstar_d_rho(rho_vec[k])*(d_pstar_d_rho(rho_vec[k]) - 1) + share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + (d_pstar_d_rho(rho_vec[k]) - 1)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]))
								term3 = (1-est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
								res = term1*int + term2*term3
								whess_mat[i,j] = -res
							else
								whess_mat[i,j] = 0.0
							end
						end
					end
					# rho lambda part of hess (upper right)
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 2:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k])) + lambda_vec[k]*(-est_pdf(lambda_vec[k]))
								res = term1 + term2*term3
								whess_mat[i,j+N-1-1] = -res
							elseif j == i+1
								term1 = m*((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*(lambda_vec[k+1]*est_pdf(lambda_vec[k+1]))
								term2 = m*((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k+1])) + lambda_vec[k+1]*(-est_pdf(lambda_vec[k+1]))
								res = term1 + term2*(-1)*term3
								whess_mat[i,j+N-1-1] = -res
									
							else
								whess_mat[i,j+N-1-1] = 0.0
							end
						end
					end
					#lambda lambda part of hess
					for i = 2:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 2:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*share(ps1))
								term2 = (-(lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k])))
								term3 = m*((rho_vec[k-1] - c)*share(ps2))
								term4 = (lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k]))
								term5 = m*((ps1 - rho_vec[k])*share(ps1) - (ps2 - rho_vec[k-1])*share(ps2))
								term6 = -est_pdf(lambda_vec[k]) + lambda_vec[k]*(-d_est_pdf(lambda_vec[k])) + (-est_pdf(lambda_vec[k]))
								res = term1*term2 + term3*term4 + term5*term6
								whess_mat[i+N-1-1,j+N-1-1] = -res
							else
								whess_mat[i+N-1-1,j+N-1-1] = 0.0
							end
						end
					end
					# lambda rho part (lower left). should be same as rho_lambda part.
					for i = 2:N-1
						for j = 1:N-1
							whess_mat[i+N-1-1,j] = whess_mat[j,i+N-1-1]
						end
					end
				end
			else
				function w_profit(theta) #theta should contain rho_1, ..., rho_n-1 ; lambda_1, lambda_2, ... , lambda_n-1
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]]
					profit = 0.0
					for i = 1:N-1
						k = i+1 # dealing with indexing
						f(l) =  l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						inc = ((rho_vec[k] - c)*m*share(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*share(ps1)*m - (ps2 - rho_vec[k-1])*share(ps2)*m)
						profit = profit + inc
					end
					return -profit
				end
				function wfocs!(theta, wfocs_vec)
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					# Calculating FOCs
					for i in 1:N-1 
						k = i+1 # dealing with the increment in indexing and making code look more like notation
						#calculating integral for rho FOC
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						#rho FOC
							term1 = ((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*int*m
							term2 = ((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))*m
							term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
							res = term1 + term2*term3
							wfocs_vec[i] = -res
						# lambda FOC
						if i == 1
							term1 = ((rho_vec[k] - c)*m*share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((ps1 - rho_vec[k])*m*share(ps1) - (ps2 - rho_vec[k-1])*m*share(ps2))
							term3 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2*term3
							wfocs_vec[i+N-1] = -res
						else
							term1 = ((rho_vec[k] - c)*m*share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
							term2 = ((rho_vec[k-1] - c)*m*share(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
							term3 = ((ps1 - rho_vec[k])*m*share(ps1) - (ps2 - rho_vec[k-1])*m*share(ps2))
							term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
							res = term1 + term2 + term3*term4
							wfocs_vec[i+N-1] = -res
						end
						
					end
				end
				function whess!(theta, whess_mat)
					lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
					rho_vec = [rho_0 ; theta[1:N-1]] # that inital value is just to force share(p_star(rho_0,lambda_0)) = 0. Ideally should be Inf
					
					#rho rho part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						f(l) = l*est_pdf(l)
						int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*(d_share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + d_pstar_d_rho(rho_vec[k])^2*dd_share(ps1)) + 2*d_share(ps1)*d_pstar_d_rho(rho_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*(d_share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + d_pstar_d_rho(rho_vec[k])^2*dd_share(ps1)) + d_share(ps1)*d_pstar_d_rho(rho_vec[k])*(d_pstar_d_rho(rho_vec[k]) - 1) + share(ps1)*d2_pstar_d2_rho(rho_vec[k]) + (d_pstar_d_rho(rho_vec[k]) - 1)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]))
								term3 = (1-est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
								res = term1*int + term2*term3
								whess_mat[i,j] = -res
							else
								whess_mat[i,j] = 0.0
							end
						end
					end
					# rho lambda part of hess (upper right)
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
								term2 = m*((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k])) + lambda_vec[k]*(-est_pdf(lambda_vec[k]))
								res = term1 + term2*term3
								whess_mat[i,j+N-1] = -res
							elseif j == i+1
								term1 = m*((rho_vec[k] - c)*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1))*(lambda_vec[k+1]*est_pdf(lambda_vec[k+1]))
								term2 = m*((ps1 - rho_vec[k])*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))
								term3 = (1-est_cdf(lambda_vec[k+1])) + lambda_vec[k+1]*(-est_pdf(lambda_vec[k+1]))
								res = term1 + term2*(-1)*term3
								whess_mat[i,j+N-1] = -res
									
							else
								whess_mat[i,j+N-1] = 0.0
							end
						end
					end
					#lambda lambda part of hess
					for i = 1:N-1
						k = i+1 #dealing with 1 indexing
						# Pre-calculating some stuff to avoid repeated calls to p_star
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
						for j = 1:N-1
							if j == i
								term1 = m*((rho_vec[k] - c)*share(ps1))
								term2 = (-(lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k])))
								term3 = m*((rho_vec[k-1] - c)*share(ps2))
								term4 = (lambda_vec[k]*d_est_pdf(lambda_vec[k]) + est_pdf(lambda_vec[k]))
								term5 = m*((ps1 - rho_vec[k])*share(ps1) - (ps2 - rho_vec[k-1])*share(ps2))
								term6 = -est_pdf(lambda_vec[k]) + lambda_vec[k]*(-d_est_pdf(lambda_vec[k])) + (-est_pdf(lambda_vec[k]))
								res = term1*term2 + term3*term4 + term5*term6
								whess_mat[i+N-1,j+N-1] = -res
							else
								whess_mat[i+N-1,j+N-1] = 0.0
							end
						end
					end
					# lambda rho part (lower left). should be same as rho_lambda part.
					for i = 1:N-1
						for j = 1:N-1
							whess_mat[i+N-1,j] = whess_mat[j,i+N-1]
						end
					end
				end
			end
			
			# checking hessian and gradient
			#=		
			innerx0 = [20.0, 10.0, 5.0, 3.0, 0.3, 0.6, 0.9]
			println(innerx0)
			gtest1 = ones(length(innerx0))
			gtest2 = ones(length(innerx0))
			htest = ones(length(innerx0),length(innerx0))
			eps = zeros(length(innerx0)) 
			step = 1e-9
			eps[2] = step
			upx = innerx0 + eps
			downx = innerx0 - eps
			est_grad = (w_profit(upx) - w_profit(downx))/(2*step)
			println("Numerical grad: ", est_grad)
			wfocs!(innerx0,gtest1)
			println(gtest1)
			
			wfocs!(upx,gtest1)
			wfocs!(downx,gtest2)
			whess!(innerx0,htest)
			println("Numerical Hessian: ", (gtest1 - gtest2)/(2*step))
			println(htest)
			=#		

			# using liner approx as hot start
			hs_sched = Lprice_sched_calc([0.0,0.0,0.0],N)
			hs_rhos = hs_sched[1]
			hs_ff = hs_sched[2]
			hs_lambdas = hs_sched[3]
			c = 0.0
			a = 1.0 # note not exp here
			b = 1.0 # note not exp here
			innerx0 = [hs_rhos[2:end] ; hs_lambdas[3:end-1]]
			res = Optim.optimize(w_profit,wfocs!,whess!,innerx0,method=NewtonTrustRegion(), extended_trace=false)
			innerx0=Optim.minimizer(res)
			steps = round(Int,ceil((params[1] - 1.0 + 1.0)/0.1)) # divide by desired step size. Could be made smaller for smoother continuation
			c_space = linspace(1.0,params[1],steps)
			a_space = linspace(1.0,exp(params[2]),steps)
			b_space = linspace(1.0,exp(params[3]),steps)
			for i = 1:steps	
				#println([c;innerx0])
				c = c_space[i]
				a = a_space[i]
				b = b_space[i]
				res = Optim.optimize(w_profit,wfocs!,whess!,innerx0,method=NewtonTrustRegion(), extended_trace=false)
				innerx0=Optim.minimizer(res)
			end
			sol_sched = innerx0	
			est_rhos = [rho_0 ; sol_sched[1:N-1]]
			est_lambdas = [lambda_lb ; 0.0 ; sol_sched[N:end] ; lambda_ub]

			
			# Calculating fixed fees
			est_ff = [0.0]
			for i in 1:N-1
				k = i+1
				A = est_ff[k-1] + (p_star(est_rhos[k]) - est_rhos[k])*share(p_star(est_rhos[k]))*M*est_lambdas[k] - (p_star(est_rhos[k-1]) - est_rhos[k-1])*share(p_star(est_rhos[k-1]))*M*est_lambdas[k]
				push!(est_ff,A)
			end
			return (est_rhos,est_ff,est_lambdas)
			
		end
		

		obs_N = length(obs_rhos)+1
		
		# testing recovery of params with fake data
		#=	
		println("Testing recovery of parameters with 'fake' data")
		x0 = [5.0; 1.0; 2.0]
		x0 = [18.0;0.0;4.0]
		nlrho,nlff,nllamb = price_sched_calc(x0,obs_N) # repeated calls give slightly different answers. Can't track down the bug
		obs_rhos = nlrho[2:end]
		obs_ff = nlff[2:end]
		obs_lambdas = nllamb[2:end-1]
		obs_sched = [obs_rhos ; obs_lambdas]
		println(obs_sched)
		=#	
		
		# Generating Deviations. Same for all params, so only do once.
		dev_step = .05
		dev_steps = [1.0+dev_step, 1.0, 1.0-dev_step]
		cart_prod_arg = fill(dev_steps,2*(obs_N-1))
		dev_cart_prod = collect(Iterators.product(cart_prod_arg...)) # need to specify Iterators. here because I use product as another symbol
		dev_cart_array = map(collect,dev_cart_prod)
		filter!(x->x!=ones(2*(obs_N-1)),dev_cart_array) # removing no-deviation vector. Just a by product of above proceedure
		dev_sched_vec =  map(x->obs_sched.*x,dev_cart_array)
		dev_pstar_vec = []
		for s in dev_sched_vec
			push!(dev_pstar_vec,map(p_star,[2*max_rho;s[1:obs_N-1]]))
		end
			
		function moment_obj_func(omega::Array{Float64,1})
			c = omega[1] # marginal cost for wholesaler
			
			a = exp(omega[2]) #1.0 # first dist param
			b = exp(omega[3]) # second dist param 
			
			lambda_lb = 0.0
			lambda_ub = 1.0
				
			est_cdf(x::Float64) = cdf(Beta(a,b),x)
			est_pdf(x::Float64) = pdf(Beta(a,b),x)
			
			rho_0 = 2.0*max_rho
			
			function w_profit(sched::Array{Real,1}, p_star_pre = Float64[])
				theta = collect(sched)
				lambda_vec = [lambda_lb; theta[obs_N:end]; lambda_ub] 
				rho_vec = [rho_0 ; theta[1:obs_N-1]]
				profit = 0.0
				for i = 1:obs_N-1
					k = i+1 # dealing with indexing
					f(l::Float64) =  l*est_pdf(l)
					int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
					# Pre-calculating some stuff to avoid repeated calls to p_star
					if length(p_star_pre) == 0
						ps1 = p_star(rho_vec[k])
						ps2 = p_star(rho_vec[k-1])
					else
						ps1 = p_star_pre[k]
						ps2 = p_star_pre[k-1] # note that p_star_re should include p_star(rho_0)
					end
					inc = ((rho_vec[k] - c)*M*share(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*share(ps1)*M - (ps2 - rho_vec[k-1])*share(ps2)*M)
					profit = profit + inc
				end
				return profit
			end
			dev_profit = map(w_profit,dev_sched_vec,dev_pstar_vec)
			obs_profit = w_profit(obs_sched)

			obs_vec = fill(obs_profit,length(dev_sched_vec))
			profit_diff = obs_vec - dev_profit
			profit_diff = convert(Array{Float64,1},profit_diff) # need to convert. Not sure why it's not typed above
			profit_diff[obs_vec .< 0.0] = -99999.0 # if profit < 0 with observed sched, then it must violate model. setting to large value
			min2_vec = min(0.0,profit_diff).^2.0
			res = sum(min2_vec)/length(dev_sched_vec)
			if (b >= exp(4.5)) | (c <= 0.0) | (a <= 1.0)
				res = Inf
			end
			return res
		end
		
		function ∇moment_obj_func(grad,params)
			eps = 1e-9
			grad[1] = (moment_obj_func(params[1]+eps,params[2],params[3]) - moment_obj_func(params[1]-eps,params[2],params[3]))/(2*eps)
			grad[2] = (moment_obj_func(params[1],params[2]+eps,params[3]) - moment_obj_func(params[1],params[2]-eps,params[3]))/(2*eps)
			grad[3] = (moment_obj_func(params[1],params[2],params[3]+eps) - moment_obj_func(params[1],params[2],params[3]-eps))/(2*eps)
		end
			
		# grid search approach	
		min_rho = minimum(obs_rhos)
		c_grid = floor(.5*min_rho):.1:ceil(1.5*min_rho)
		b_grid = -4:.05:4
		println(c_grid)
		println(b_grid)
		cart_grid = Iterators.product(c_grid,b_grid)
		#=
		Q_eval = Array{Float64}(0,3)
		for c in c_grid
			for b in b_grid
				q = moment_obj_func(c,b)
				res = [c b q]
				Q_eval = [Q_eval ; res]
			end
		end
		=#
		#Q_eval = map((args)->moment_obj_func(args...),cart_grid) # grid search - potentially able to parallelize this with pmap
		

		# Simulated annealing approach
		sa_temp(t) = 0.0
		hs_res = Optim.optimize(moment_obj_func, [1.0,1.0,1.0],method=NelderMead(parameters= Optim.FixedParameters()), iterations = 10000, extended_trace = true, store_trace = true, show_every = false)
		println(hs_res)
		x0 = Optim.minimizer(hs_res)
		res = Optim.optimize(moment_obj_func, x0,method = SimulatedAnnealing(temperature=sa_temp), iterations = 10000, extended_trace = true, store_trace = true, show_every = false)
		println(res)
		Q_eval = Optim.f_trace(res)
		x_eval = Optim.x_trace(res)

		csvfile = open("moment_$market"*"_$product.csv", "w")
		write(csvfile, "c,a,b,Q\n")
		out_mat = [[y[i] for y in x_eval, i in 1:3] Q_eval]
		writedlm(csvfile,out_mat,',')
		close(csvfile)
		
		
		#min_q_ind = findin(Q_eval[:,3],minimum(Q_eval[:,3]))
		min_q = minimum(Q_eval)
		min_q_ind = find(Q_eval .<= min_q+1e-4)
		feas_c = out_mat[min_q_ind,1]
		feas_a = out_mat[min_q_ind,2]
		feas_b = out_mat[min_q_ind,3]
		c_ub = maximum(feas_c)
		c_lb = minimum(feas_c)
		a_ub = maximum(feas_a)
		a_lb = minimum(feas_a)
		b_ub = maximum(feas_b)
		b_lb = minimum(feas_b)
		println([c_lb,c_ub])
		println([a_lb,a_ub])
		println([b_lb,b_ub])
			
		
		# finding xi param (cost of additional segment)
		println(M)	
		mid_c = (c_lb + c_ub)/2.0
		mid_a = (a_lb + a_ub)/2.0
		mid_b = (b_lb + b_ub)/2.0
		println([mid_c,mid_b])
		eq_ps = price_sched_calc([mid_c,mid_a,mid_b], obs_N)
		println(eq_ps)
		more_ps = price_sched_calc([mid_c,mid_a,mid_b], obs_N+1)
		println(more_ps)
		less_ps = price_sched_calc([mid_c,mid_a,mid_b],obs_N-1)
		println(less_ps)
		
		println(M)	
		less_rho = less_ps[1]
		less_ff = less_ps[2]
		less_lambda = less_ps[3]
		less_sched = [less_rho[2:end] ; less_lambda[2:end-1]]
		more_rho = more_ps[1]
		more_ff = more_ps[2]
		more_lambda = more_ps[3]
		more_sched = [more_rho[2:end] ; more_lambda[2:end-1]]
		eq_rho = eq_ps[1]
		eq_ff = eq_ps[2]
		eq_lambda = eq_ps[3]
		eq_sched = [eq_rho[2:end] ; eq_lambda[2:end-1]]
		
		function out_w_profit(sched, p_star_pre = Float64[])
			c = mid_c # marginal cost for wholesaler
			
			a = exp(mid_a) # first dist param
			b = exp(mid_b) # second dist param 
			
			lambda_lb = 0.0
			lambda_ub = 1.0

			N = div(length(sched),2) + 1 # integer division
				
			est_cdf(x) = cdf(Beta(a,b),x)
			est_pdf(x) = pdf(Beta(a,b),x)
			
			rho_0 = 2.0*max_rho
			theta = collect(sched)
			lambda_vec = [lambda_lb; theta[N:end]; lambda_ub] 
			rho_vec = [rho_0 ; theta[1:N-1]]
			profit = 0.0
			for i = 1:N-1
				k = i+1 # dealing with indexing
				f(l) =  l*est_pdf(l)
				int = sparse_int(f,lambda_vec[k],lambda_vec[k+1])
				# Pre-calculating some stuff to avoid repeated calls to p_star
				if length(p_star_pre) == 0
					ps1 = p_star(rho_vec[k])
					ps2 = p_star(rho_vec[k-1])
				else
					ps1 = p_star_pre[k]
					ps2 = p_star_pre[k-1] # note that p_star_re should include p_star(rho_0)
				end
				inc = ((rho_vec[k] - c)*M*share(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*share(ps1)*M - (ps2 - rho_vec[k-1])*share(ps2)*M)
				profit = profit + inc
			end
			return profit
		end
		eq_wprofit = out_w_profit(eq_sched)
		less_wprofit = out_w_profit(less_sched)
		more_wprofit = out_w_profit(more_sched)
		xi_ub = eq_wprofit - less_wprofit
		xi_lb = more_wprofit - eq_wprofit
		println([xi_lb,xi_ub])
		mid_xi = (xi_lb + xi_ub)/2
		println(mid_xi)
		
		# Linear price counter factual
		lin_ps = price_sched_calc([mid_c,mid_a,mid_b],2)
		println(lin_ps)
		lin_rho = lin_ps[1]
		lin_ff = lin_ps[2]
		lin_lambda = lin_ps[3]
		lin_sched = [lin_rho[2:end] ; lin_lambda[2:end-1]]
		lin_wprofit = out_w_profit(lin_sched)
	
		function rprofit(rho,lambda,ff)
			ps = p_star(rho)
			return (ps - rho)*share(ps)*M*lambda - ff
		end
		eq_rprofit = 0.0
		for i = 1:obs_N-1
			k = i + 1
			f(x) = rprofit(eq_rho[k],x,eq_ff[k])
			res = sparse_int(f,eq_lambda[k],eq_lambda[k+1])
			eq_rprofit += res
		end

		eq_csurp = 0.0
		tot_q = 0.0
		for i = 1:obs_N-1
			k = i + 1
			ps = p_star(eq_rho[k])
			q = share(ps)
			res = quadgk(share,ps,Inf)[1]
			eq_csurp += res*q
			tot_q += q
		end
		eq_csurp = eq_csurp / tot_q

		lin_ps = p_star(lin_rho[2])
		println("Linear retail price: ", lin_ps)
		lin_csurp = quadgk(share,lin_ps,Inf)[1]


		lrp(x) = rprofit(lin_rho[2],x,lin_ff[2])
		lin_rprofit = sparse_int(lrp,lin_lambda[2],lin_lambda[3])
	
			
		Δw_profit = ((lin_wprofit-2*mid_xi) - (eq_wprofit - obs_N*mid_xi))/(eq_wprofit - obs_N*mid_xi)
		Δr_profit = (lin_rprofit - eq_rprofit)/(eq_rprofit)
		Δc_surp = (lin_csurp - eq_csurp)/(eq_csurp)
		println("Change in wholesaler profit: ",Δw_profit)
		println("Change in retailer profit: ",Δr_profit)
		println("Change in consumer surplus: ", Δc_surp)
		
		return 1
	else
		println("Product has no matching price data.")
	end
end

df = readtable("../../demand_estimation/berry_logit/berry_logit.csv")

tups = [] # market,product tuples to run estimation on

#markets = convert(Vector, levels(df[:,:mkt]))
markets = [11725]
for market in markets
	#products = levels(df[df[:mkt] .== market, :product])
	products = [350]
	for product in products
		mkt = convert(Int,market)
		pd = convert(Int,product)
		push!(tups,(mkt,pd))
	end
end

total_tups = length(tups)
println("Total estimates to produce: $total_tups")
ss_est((11725,350))
#para_out = pmap(ss_est,tups)
