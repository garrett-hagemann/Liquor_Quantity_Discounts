@everywhere using DataArrays, DataFrames, ForwardDiff, NLsolve, Roots, Distributions, Optim, NLopt, JuMP

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

@everywhere function sparse_int(f::Function, a::Number, b::Number)
	#= Implements sparse grid quadrature from sparsegrids.de
	This implements the 1 dimensional rule that is exact for
	polynomials up to order 25.

	f = one dimensional function
	a = lower bound of integration
	b = upperbound of integration

	=#
	 weights = [.052328105232810528, .13424401342440137, .20069872006987202, .22545832254583226, .20069872006987202, .13424401342440137, .052328105232810528]

	nodes = [.01975439999999995,.11270170000000002,.28287810000000002,.5,.71712189999999998, .88729829999999998,.98024560000000005]

	f_evals = [f((b-a)*u + a) for u in nodes]
	return dot(f_evals, weights)*(b-a)::Float64
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
		#M = df[prod_bool,:M][1]
		M = 1.0
		prod_price = df[prod_bool, :price][1]

		# lists of observed prices and quantities. Need these to calculate error using 
		# estimated price schedule 
		rho_list = [:actual_p0, :actual_p1, :actual_p2, :actual_p3, :actual_p4, :actual_p5, :actual_p6,
			:actual_p7, :actual_p8, :actual_p9, :actual_p10]

		cutoff_q_list = [:disc_q0, :disc_q1, :disc_q2, :disc_q3, :disc_q4, :disc_q5, :disc_q6,
			:disc_q7, :disc_q8, :disc_q9, :disc_q10]

		obs_rhos = dropna(convert(DataArray,df[prod_bool, rho_list])'[:,1]) #some goofy conversion to make the dropna work as expected
		obs_cutoff_q = dropna(convert(DataArray,df[prod_bool, cutoff_q_list])'[:,1]) #same as above
		obs_ff = [0.0] # first fixed fee is by definition
		
		max_rho = copy(maximum(obs_rhos))
		max_mc = 15.0
		# calculating series of A (intercepts of piecewise linear price schedules
		# calculating tariff
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

		### NLsolve pstar
		#=	
		function p_star(rho,l)
			function g!(p,gvec)
				gvec[1] =  (p - rho + l)*d_share(p)*M + share(p)*M
			end
			function gj!(p, gjvec)
				gjvec[1] = (p - rho + l)*dd_share(p)*M + 2.0*d_share(p)*M
			end
			res = nlsolve(g!,gj!,[rho-l],show_trace = false, extended_trace = false, method = :trust_region) 
			return res.zero[1]
		end
		=#	
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

			# params are the params of the wholesaler's problem we're trying to
			# estimate. 

			# Function returns the coefficients defining the price schedule as [rho, lambda]
			# where rho is n long and lambda is n-1 long for an n option schedule

			c = params[1] # marginal cost for wholesaler
			#max_mc = params[2] # max MC for retailer. Scales type parameter
			#M = exp(params[5])
			#distribution parameters
			a = 1.0 # exp(params[2]) #
			b = exp(params[2]) # 
			lambda_lb = 0.0
			lambda_ub = 1.0
				
			est_cdf(x) = cdf(Beta(a,b),x)
			est_pdf(x) = pdf(Beta(a,b),x)
		
			if (constrained == 1)	
				function urho(k::Integer)
					return -LA/LB + (LA + LB*c)/(2*LB)*(2*k - 1)/(N-1)
				end
				function ulambda(k::Integer)
					return (k-1)/(N-1)
				end
				rho_0 = urho(0)
			end
			##### Optimizer Approach
			function Lw_profit(params...)
				theta = collect(params)
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
					inc = ((rho_vec[k] - c)*M*Lshare(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*Lshare(ps1)*M - (ps2 - rho_vec[k-1])*Lshare(ps2)*M)
					profit = profit + inc
				end
				return -profit
			end
			function Lwfocs!(wfocs_vec,params...)
				theta = collect(params)
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
						term1 = ((rho_vec[k] - c)*M*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + M*Lshare(ps1))*int
						term2 = ((ps1 - rho_vec[k])*M*Ld_share(ps1)*Ld_pstar_d_rho(rho_vec[k]) + M*Lshare(ps1)*(Ld_pstar_d_rho(rho_vec[k]) - 1))
						term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
						res = term1 + term2*term3
						wfocs_vec[i] = -res
					# lambda FOC
					if i == 1
						term1 = ((rho_vec[k] - c)*M*Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
						term2 = ((ps1 - rho_vec[k])*M*Lshare(ps1) - (ps2 - rho_vec[k-1])*M*Lshare(ps2))
						term3 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
						res = term1 + term2*term3
						wfocs_vec[i+N-1] = -res
					else
						term1 = ((rho_vec[k] - c)*M*Lshare(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
						term2 = ((rho_vec[k-1] - c)*M*Lshare(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
						term3 = ((ps1 - rho_vec[k])*M*Lshare(ps1) - (ps2 - rho_vec[k-1])*M*Lshare(ps2))
						term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
						res = term1 + term2 + term3*term4
						wfocs_vec[i+N-1] = -res
					end
					
				end
			end
			rho_start = convert(Array{Float64,1},[urho(k) for k = 1:N-1])
			lambda_start = convert(Array{Float64,1},[ulambda(k) for k = 1:N-1])	
			innerx0 = [rho_start;lambda_start]
			#innerx0 = [1.0, 0.75, 0.5, 0.3, 0.6, 0.9] + 0*randn(2*(N-1))
			# checking hessian and gradient
			#=
			println(innerx0)
			
			gtest1 = ones(2*(N-1))
			gtest2 = ones(2*(N-1))
			htest = ones(2*(N-1),2*(N-1))
			eps = zeros(2*(N-1)) 
			step = 1e-9
			eps[4] = step
			est_grad = (Lw_profit(innerx0+eps) - Lw_profit(innerx0-eps))/(2*step)
			println("Numerical grad: ", est_grad)
			Lwfocs!(innerx0,gtest1)
			println(gtest1)

			solution = 1
			return solution
			=#
			# Solving problem numerically. Using JuMP modeling language
			try
				JuMP.register(:Lshare,1,Lshare,autodiff = true) # share function
			catch e
				if e.msg == "Operator Lshare has already been defined"
					ind = pop!(ReverseDiffSparse.univariate_operator_to_id, :Lshare)
					deleteat!(ReverseDiffSparse.univariate_operators,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_f,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprime,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprimeprime,ind)
					JuMP.register(:Lshare,1,Lshare,autodiff = true) # share function
				end
			end

			try
				JuMP.register(:Lp_star,1,Lp_star,Ld_pstar_d_rho,Ld2_pstar_d2_rho) 
			catch e
				if e.msg == "Operator Lp_star has already been defined"
					ind = pop!(ReverseDiffSparse.univariate_operator_to_id, :Lp_star)
					deleteat!(ReverseDiffSparse.univariate_operators,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_f,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprime,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprimeprime,ind)
					JuMP.register(:Lp_star,1,Lp_star,Ld_pstar_d_rho,Ld2_pstar_d2_rho) 
				end
			end

			try
				JuMP.register(:Lw_profit,2*(N-1),Lw_profit,Lwfocs!, autodiff = false)
			catch e
				if e.msg == "Operator Lw_profit has already been defined"
					ind = pop!(ReverseDiffSparse.operator_to_id, :Lw_profit)
					pop!(ReverseDiffSparse.user_operator_map,ind)
					JuMP.register(:Lw_profit,2*(N-1),Lw_profit,Lwfocs!, autodiff = false)
				end
			end
			global Linner_m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP, ftol_abs=1e-6, ftol_rel=1e-6, maxeval = 1000)) #bad form, but needed for meta programming BS
			@variable(Linner_m,Linner_s[1:2*(N-1)])
			global Linner_s = Linner_s #bad form, but needed for meta programming BS
			for k = 1:N-1 
				setvalue(Linner_s[k],urho(k)) # setting rho starting values
				setvalue(Linner_s[N-1+k],ulambda(k)) # setting lambda starting values
			end
			if constrained == 1
				@NLconstraint(Linner_m,cons1,Linner_s[N]*Lshare(Lp_star(Linner_s[1])) == 0)
			end
		
			# defining objective function. Need to use meta-programming BS to get around syntax limitations
			obj_str = "@NLobjective(Linner_m, Min, Lw_profit("
			for i = 1:2*(N-1)
				obj_str = obj_str * "Linner_s[$i],"
			end
			obj_str = obj_str[1:end-1] # removing last comma
			obj_str = obj_str * "))" # closing parens
			eval(parse(obj_str))
			inner_sol = solve(Linner_m)
			sol_sched = getvalue(Linner_s)
			
			est_rhos = [rho_0 ; sol_sched[1:N-1]]
			est_lambdas = [lambda_lb ; sol_sched[N:end] ; lambda_ub]

			
			# Calculating fixed fees
			est_ff = [0.0]
			for i in 1:N-1
				k = i+1
				A = est_ff[k-1] + (Lp_star(est_rhos[k]) - est_rhos[k])*Lshare(Lp_star(est_rhos[k]))*M*est_lambdas[k] - (Lp_star(est_rhos[k-1]) - est_rhos[k-1])*Lshare(Lp_star(est_rhos[k-1]))*M*est_lambdas[k]
				push!(est_ff,A)
			end
			return (est_rhos,est_ff,est_lambdas)
		end		
		
		function price_sched_calc(params,N)
			constrained = 1 # solving problem with A1 constrained to 0

			# params are the params of the wholesaler's problem we're trying to
			# estimate. 

			# Function returns the coefficients defining the price schedule as [rho, lambda]
			# where rho is n long and lambda is n-1 long for an n option schedule

			c = params[1] # marginal cost for wholesaler
			#max_mc = params[2] # max MC for retailer. Scales type parameter
			#M = exp(params[5])
			#distribution parameters
			a = 1.0 # exp(params[2]) #
			b = exp(params[2]) # 
			lambda_lb = 0.0
			lambda_ub = 1.0
				
			est_cdf(x) = cdf(Beta(a,b),x)
			est_pdf(x) = pdf(Beta(a,b),x)
			
			rho_0 = 2.0*max_rho
		
			##### Optimizer Approach
			function w_profit(params...)
				theta = collect(params)
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
					inc = ((rho_vec[k] - c)*M*share(ps1))*int + (1-est_cdf(lambda_vec[k]))*lambda_vec[k]*((ps1 - rho_vec[k])*share(ps1)*M - (ps2 - rho_vec[k-1])*share(ps2)*M)
					profit = profit + inc
				end
				return -profit
			end
			function wfocs!(wfocs_vec,params...)
				theta = collect(params)
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
						term1 = ((rho_vec[k] - c)*M*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + M*share(ps1))*int
						term2 = ((ps1 - rho_vec[k])*M*d_share(ps1)*d_pstar_d_rho(rho_vec[k]) + M*share(ps1)*(d_pstar_d_rho(rho_vec[k]) - 1))
						term3 = (1 - est_cdf(lambda_vec[k]))*lambda_vec[k] - (1-est_cdf(lambda_vec[k+1]))*lambda_vec[k+1]
						res = term1 + term2*term3
						wfocs_vec[i] = -res
					# lambda FOC
					if i == 1
						term1 = ((rho_vec[k] - c)*M*share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
						term2 = ((ps1 - rho_vec[k])*M*share(ps1) - (ps2 - rho_vec[k-1])*M*share(ps2))
						term3 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
						res = term1 + term2*term3
						wfocs_vec[i+N-1] = -res
					else
						term1 = ((rho_vec[k] - c)*M*share(ps1))*(-lambda_vec[k]*est_pdf(lambda_vec[k]))
						term2 = ((rho_vec[k-1] - c)*M*share(ps2))*(lambda_vec[k]*est_pdf(lambda_vec[k]))
						term3 = ((ps1 - rho_vec[k])*M*share(ps1) - (ps2 - rho_vec[k-1])*M*share(ps2))
						term4 = (1 - est_cdf(lambda_vec[k]) - lambda_vec[k]*est_pdf(lambda_vec[k]))
						res = term1 + term2 + term3*term4
						wfocs_vec[i+N-1] = -res
					end
					
				end
			end
			innerx0 = [1.0, 0.75, 0.5, 0.3, 0.6, 0.9]
			# checking hessian and gradient
			#=		
			println(innerx0)
			
			gtest1 = ones(2*(N-1))
			gtest2 = ones(2*(N-1))
			htest = ones(2*(N-1),2*(N-1))
			eps = zeros(2*(N-1)) 
			step = 1e-9
			eps[6] = step
			upx = innerx0 + eps
			downx = innerx0 - eps
			est_grad = (w_profit(upx...) - w_profit(downx...))/(2*step)
			println("Numerical grad: ", est_grad)
			wfocs!(gtest1,innerx0...)
			println(gtest1)

			solution = 1
			return solution
			=#
			# Solving problem numerically. Using JuMP modeling language
			try
				JuMP.register(:share,1,share,autodiff = true) # share function
			catch e
				if e.msg == "Operator share has already been defined"
					ind = pop!(ReverseDiffSparse.univariate_operator_to_id, :share)
					deleteat!(ReverseDiffSparse.univariate_operators,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_f,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprime,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprimeprime,ind)
					JuMP.register(:share,1,share,autodiff = true) # share function
				end
			end

			try
				JuMP.register(:p_star,1,p_star,d_pstar_d_rho,d2_pstar_d2_rho) 
			catch e
				if e.msg == "Operator p_star has already been defined"
					ind = pop!(ReverseDiffSparse.univariate_operator_to_id, :p_star)
					deleteat!(ReverseDiffSparse.univariate_operators,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_f,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprime,ind)
					pop!(ReverseDiffSparse.user_univariate_operator_fprimeprime,ind)
					JuMP.register(:p_star,1,p_star,d_pstar_d_rho,d2_pstar_d2_rho) 
				end
			end

			try
				JuMP.register(:w_profit,2*(N-1),w_profit,wfocs!, autodiff = false)
			catch e
				if e.msg == "Operator w_profit has already been defined"
					ind = pop!(ReverseDiffSparse.operator_to_id, :w_profit)
					pop!(ReverseDiffSparse.user_operator_map,ind)
					JuMP.register(:w_profit,2*(N-1),w_profit,wfocs!, autodiff = false)
				end
			end


			global inner_m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP, ftol_abs=1e-6, ftol_rel=1e-6, maxeval=1000)) #bad form, but needed for meta programming BS
			@variable(inner_m,inner_s[1:2*(N-1)])
			global inner_s = inner_s #bad form, but needed for meta programming BS
			
			# setting starting values. Using solution to linear model as hot start
			(Lrho,Lff,Llambda) = Lprice_sched_calc(params,N)
			for k = 1:N-1 
				setvalue(inner_s[k],Lrho[k+1]) # setting rho starting values
				setvalue(inner_s[N-1+k],Llambda[k+1]) # setting lambda starting values
			end
			if constrained == 1
				@NLconstraint(inner_m,cons1,inner_s[N]*share(p_star(inner_s[1])) == 0)
			end
		
			# defining objective function. Need to use meta-programming BS to get around syntax limitations
			obj_str = "@NLobjective(inner_m, Min, w_profit("
			for i = 1:2*(N-1)
				obj_str = obj_str * "inner_s[$i],"
			end
			obj_str = obj_str[1:end-1] # removing last comma
			obj_str = obj_str * "))" # closing parens
			eval(parse(obj_str))
			inner_sol = solve(inner_m)
			sol_sched = getvalue(inner_s)
			
			est_rhos = [rho_0 ; sol_sched[1:N-1]]
			est_lambdas = [lambda_lb ; sol_sched[N:end] ; lambda_ub]

			
			# Calculating fixed fees
			est_ff = [0.0]
			for i in 1:N-1
				k = i+1
				A = est_ff[k-1] + (p_star(est_rhos[k]) - est_rhos[k])*share(p_star(est_rhos[k]))*M*est_lambdas[k] - (p_star(est_rhos[k-1]) - est_rhos[k-1])*share(p_star(est_rhos[k-1]))*M*est_lambdas[k]
				push!(est_ff,A)
			end
			return (est_rhos,est_ff,est_lambdas)
			
		end
		function obj_func(omega::Vector, N::Int, W::Matrix)
			rho_hat,ff_hat,lambda_hat = price_sched_calc(omega,N)
			vec = [(rho_hat[2:end] - obs_rhos) ; (ff_hat[3:end] - obs_ff[2:end])]'
			res = vec*W*vec'
			return res[1]
		end

		N = length(obs_rhos)+1
		#W = eye(2*(N-1)- 1)
		# testing recovery of params with fake data
		#=
		println("Testing recovery of parameters with 'fake' data")
		x0 = [15.0; log(5.0)]
		nlrho,nlff,nllamb = price_sched_calc(x0,N)
		obs_rhos = nlrho[2:end]
		obs_ff = [0.0;nlff[3:end]]
		println(obs_rhos)
		=#


		#ux0 = [3.0,15.0]
		W = Diagonal([1./(obs_rhos.^2) ; 1./(obs_ff[2:end].^2)])*eye(2*N-3)
		# checking Objective func gradient
		#eps = zeros(2)
		#eps[1] = 1e-9
		#println((uLobj_func(ux0+eps,N,W) - uLobj_func(ux0-eps,N,W))/(2*1e-9))
		#jtest = zeros(2,1)
		#uLobj_foc!(ux0,jtest,N,W)	
		#println(jtest)
		
		# testing hot start
		#x0 = [0.0; 2.0; log(1.0); log(1.0)]
		#hsrho,hsff,hslambda = Lprice_sched_calc(x0,N)
		#hs = [hsrho[2:end],hslambda[2:end-1]]
		#println(Lprice_sched_calc(x0,5))
		#println(price_sched_calc(x0,4))
		#println(price_sched_calc(x0,N; hot_start = hs))
		
		# Optimizing uniform linear model	
		#uLg(x) = uLobj_func(x,N,W)
		#focs!(x,vec) = uLobj_foc!(x,vec,N,W)
		#optim_res = optimize(uLg,focs!,ux0,BFGS(),OptimizationOptions(show_every = false, extended_trace = false, iterations = 1500, g_tol = 1e-6))
		#println(optim_res)
		#min_X = Optim.minimizer(optim_res)
		
		outerg(x) = obj_func(x,N,W)
			
		# grid search
		#grid_points = [[i,j,0.0,0.0] for i = 0.0:10.0, j = 2.0:2.0:20.0]
		println("Starting initial grid search")
		grid_points = [[i,l] for i = linspace(0,20,20), l=linspace(0.0,2.0,5)]
		grid_evals = map(outerg,grid_points)
		min_ind = indmin(grid_evals)
		x0 = grid_points[min_ind]
		println("Best Starting Guess from grid: ",x0)
		
		outerg(x) = obj_func(x,N,W)
		#x0 = [2.0; 11.0; log(1.0); log(1.0)]
		optim_res = Optim.optimize(outerg,x0,NelderMead(),OptimizationOptions(show_every = false, extended_trace = false, iterations = 1500, g_tol = 1e-8))
		println(optim_res)
		min_X = Optim.minimizer(optim_res)
		hsrho,hsff,hslambda = Lprice_sched_calc(min_X,N)
		hs = [hsrho[2:end];hslambda[2:end-1]]
		fit_ps = price_sched_calc(min_X,N)

		outtuple = (product,market,min_X[1], max_mc, 1.0, min_X[2], fit_ps)
		return outtuple
		
	else
		println("Product has no matching price data.")
	end
end

df = readtable("../../demand_estimation/berry_logit/berry_logit.csv")

tups = [] # market,product tuples to run estimation on

markets = convert(Vector, levels(df[:,:mkt]))
#markets = [11725]
for market in markets
	products = levels(df[df[:mkt] .== market, :product])
	#products = [350]
	for product in products
		m = convert(Int,market)
		p = convert(Int,product)
		push!(tups,(m,p))
	end
end

total_tups = length(tups)
println("Total estimates to produce: $total_tups")
csvfile = open("indirect_est.csv", "w")
write(csvfile, "product,mkt,c,lambda_ub,a,b,price_sched\n")
para_out = pmap(ss_est,tups)
writedlm(csvfile,para_out,"|")
close(csvfile)
