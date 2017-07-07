using ForwardDiff, Distributions, Roots, Optim, NLsolve

function sparse_int(f::Function, a::Number, b::Number)
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
	return dot(f_evals, weights)*(b-a)
end

function d(p,t)
	return t*(1-p)
end

function dd_dt(p,t)
	return (1-p)
end

function dd_dp(p,t)
	return -t
end

function dd_dp_dp(p,t)
	return 0.0
end

function dd_dt_dp(p,t)
	return -1.0
end

function dd_dp_dt(p,t)
	return -1.0
end

function dd_dt_dt(p,t)
	return 0.0
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


a = 4.0
b = 4.0
c = 0.25


#t_pdf(x) = pdf(Beta(a,b),x)
#t_cdf(x) = cdf(Beta(a,b),x)
t_pdf(x) = ks_dens(x,a,b)
t_cdf(x) = ks_dist(x,a,b)

t_lb = 0.0
t_ub = 1.0

N = 3

function focs!(theta,focs_vec)

	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; 0.0; theta[N:end] ; t_ub]

	#price focs
	for i = 2:N # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i] - c)*dd_dp(rho_vec[i],t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(rho_vec[i],t))
		res_rho = sparse_int(f,t_vec[i],t_vec[i+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i]) + (1-t_cdf(t_vec[i]))/t_pdf(t_vec[i])*dd_dt(p,t_vec[i]))
		res_t = sparse_int(g,rho_vec[i],rho_vec[i-1])
		focs_vec[i-1] = res_rho
		focs_vec[i+N-2-1] = res_t # actually incorrect for i = 1. Correctly overwrites it on next iteration. Thus, wouldn't work when N = 2
	end
end

function focs_jac!(theta,focs_jac_mat)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; 0.0 ; theta[N:end] ; t_ub]
	#rho rho part of jac
	for i = 1:N-1 # i refers to scheduel part, not index position
		for j = 1:N-1
			if j == i
				f(t) = ((rho_vec[i+1] - c)*dd_dp_dp(rho_vec[i+1],t) + dd_dp(rho_vec[i+1],t) + (1-t_cdf(t))/(t_pdf(t))*(dd_dt_dp(rho_vec[i+1],t)))*t_pdf(t)
				res = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
				focs_jac_mat[i,j] = res
			else
				res = 0
				focs_jac_mat[i,j] = res
			end
		end
	end

	# rho t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				res = -((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1])/t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1],t_vec[i+1]))*t_pdf(t_vec[i+1])
				focs_jac_mat[i,j+N-1-1] = res
			elseif j == i+1
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1+1]) + (1 - t_cdf(t_vec[i+1+1]))/(t_pdf(t_vec[i+1+1]))*dd_dt(rho_vec[i+1],t_vec[i+1+1]))*t_pdf(t_vec[i+1+1])
				focs_jac_mat[i,j+N-1-1] = res
			else
				res = 0
				focs_jac_mat[i,j+N-1-1] = res
			end
		end
	end
	# t rho part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i - 1
				res = ((rho_vec[i+1-1] - c)*dd_dp(rho_vec[i+1-1],t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1-1],t_vec[i+1]))
				focs_jac_mat[i+N-1-1,j] = res
			elseif j == i
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(rho_vec[i+1],t_vec[i+1]))
				focs_jac_mat[i+N-1-1,j] = res
			else
				res = 0
				focs_jac_mat[i+N-1-1,j] = res
			end
		end
	end
	# t t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				g(p) = (p - c)*dd_dp_dt(p,t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt_dt(p,t_vec[i+1]) + dd_dt(p,t_vec[i+1])*(-1.0)
				res = sparse_int(g,rho_vec[i+1], rho_vec[i+1-1])
				focs_jac_mat[i+N-1-1,j+N-1-1] = res
			else
				res = 0
				focs_jac_mat[i+N-1-1,j+N-1-1] = res
			end
		end
	end
end

srand(200096868)
##### NLsolve solution.
#=
x0 = [3/4,1/4,1/2] + randn(3)*0.10 # works with .01, .1, but 1 is too far from truth
#test = similar(x0)
#focs!(x0,test)
#println(test)

#solution = nlsolve(focs!,x0,show_trace = true, extended_trace = true, iterations = 500)
@time solution = nlsolve(focs!,focs_jac!,x0,show_trace = false, extended_trace = false, iterations = 500, ftol = 1e-8, method = :trust_region)
@time solution = nlsolve(focs!,focs_jac!,x0,show_trace = false, extended_trace = false, iterations = 500, ftol = 1e-8, method = :trust_region)
@time solution = nlsolve(focs!,focs_jac!,x0,show_trace = false, extended_trace = false, iterations = 500, ftol = 1e-8, method = :trust_region)
@time solution = nlsolve(focs!,focs_jac!,x0,show_trace = false, extended_trace = false, iterations = 500, ftol = 1e-8, method = :trust_region)
println(solution)
=#
#=
### Implementing homotopy to get better starting guess
s = 0:.001:1
x_star = [0.6, 0.2, 0.2, 0.6]*0.0 + randn(4)*0.1
function F(theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	res = zeros(2*(N-1))
	#price focs
	for i = 2:N # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i] - c)*dd_dp(rho_vec[i],t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(rho_vec[i],t))
		res_rho = sparse_int(f,t_vec[i],t_vec[i+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i]) + (1-t_cdf(t_vec[i]))/t_pdf(t_vec[i])*dd_dt(p,t_vec[i]))
		res_t = sparse_int(g,rho_vec[i],rho_vec[i-1])
		res[i-1] = res_rho
		res[i+N-2] = res_t
	end
	return res
end
function homof!(theta,storage,s)
	res = F(theta) + (s-1)*F(x_star)
	for i in 1:length(storage)
		storage[i] = res[i]
	end
end
x0 = x_star
for i = 1:length(s)
	println("s = ",s[i])
	nlfunc(theta,storage) = homof!(theta,storage,s[i])
	solution = nlsolve(nlfunc,x0,method=:newton, show_trace=false)
	x0 = solution.zero
	println(x0)
end
println(x0)
solution = nlsolve(focs!,focs_jac!,x0)
println(solution)
=#
#=
##### Recursive Solution
function rec_solve()
	t_guess = [0; 0.0; randn(N-1-1); 1]
	t_prime_guess = [0.0; zeros(N-1);1.0]
	rho_guess = [1.0;ones(N-1)]
	eps = 1e-8
	diff = 1
	counter = 1
	tic()
	while diff > eps
		for i = 1:N-1
			function rho_hat(r)
				f(t) = ((r - c)*dd_dp(r,t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(r,t))
				res_rho = sparse_int(f,t_guess[i+1],t_guess[i+1+1])
				return res_rho
			end
			res =  fzero(rho_hat,0.0,1.0)
			rho_guess[i+1] = res
		end
		for i = 2:N-1
			function t_hat(t)
				g(p) = ((p - c)*dd_dp(p,t) + (1-t_cdf(t))/t_pdf(t)*dd_dt(p,t))
				res_t = sparse_int(g,rho_guess[i+1],rho_guess[i+1-1])
				return res_t
			end
			res = fzero(t_hat,0.0,1.0)
			t_prime_guess[i+1] = res
		end
		diff = sum(abs((t_guess - t_prime_guess)))
		t_guess = copy(t_prime_guess)
		counter = counter +1
	end
	toc()
	println(rho_guess)
	println(t_guess)
	println(counter)
end

@time rec_solve()
@time rec_solve()
@time rec_solve()
@time rec_solve()
=#

##### optimizer solution
# Unconstrained problem
function obj_func(theta)
	if !all((0 .<= theta .<= 1)) # if not all elements in [0,1]
		return Inf
	else
		rho_vec = [1; theta[1:N-1]]
		t_vec = [t_lb; theta[N:end] ; t_ub]
		profit = 0.0
		for i = 1:N-1
			f(t) = (rho_vec[i+1] - c)*d(rho_vec[i+1],t)
			int1 = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
			g(p) = d(p,t_vec[i+1])
			int2 = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
			inc = int1 + (1-t_cdf(t_vec[i+1]))*int2
			profit = profit + inc
		end
		return -profit
	end
end

function grad!(g_vec,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	#price focs
	for i = 1:N-1 # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(rho_vec[i+1],t))
		res_rho = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(p,t_vec[i+1]))
		res_t = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		g_vec[i] = -res_rho
		g_vec[i+N-1] = -res_t
	end
end
function hess!(hess_mat,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	#rho rho part of jac
	for i = 1:N-1 # i refers to scheduel part, not index position
		for j = 1:N-1
			if j == i
				f(t) = ((rho_vec[i+1] - c)*dd_dp_dp(rho_vec[i+1],t) + dd_dp(rho_vec[i+1],t) + (1-t_cdf(t))/(t_pdf(t))*(dd_dt_dp(rho_vec[i+1],t)))*t_pdf(t)
				res = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
				hess_mat[i,j] = -res
			else
				res = 0
				hess_mat[i,j] = -res
			end
		end
	end

	# rho t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				res = -((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1])/t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1],t_vec[i+1]))*t_pdf(t_vec[i+1])
				hess_mat[i,j+N-1] = -res
			elseif j == i+1
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1+1]) + (1 - t_cdf(t_vec[i+1+1]))/(t_pdf(t_vec[i+1+1]))*dd_dt(rho_vec[i+1],t_vec[i+1+1]))*t_pdf(t_vec[i+1+1])
				hess_mat[i,j+N-1] = -res
			else
				res = 0
				hess_mat[i,j+N-1] = -res
			end
		end
	end
	# t rho part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i - 1
				res = ((rho_vec[i+1-1] - c)*dd_dp(rho_vec[i+1-1],t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1-1],t_vec[i+1]))
				hess_mat[i+N-1,j] = -res
			elseif j == i
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(rho_vec[i+1],t_vec[i+1]))
				hess_mat[i+N-1,j] = -res
			else
				res = 0
				hess_mat[i+N-1,j] = -res
			end
		end
	end
	# t t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				g(p) = (p - c)*dd_dp_dt(p,t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt_dt(p,t_vec[i+1]) + dd_dt(p,t_vec[i+1])*(-1.0)
				res = sparse_int(g,rho_vec[i+1], rho_vec[i+1-1])
				hess_mat[i+N-1,j+N-1] = -res
			else
				res = 0
				hess_mat[i+N-1,j+N-1] = -res
			end
		end
	end
end
x0 = rand(2*(N-1))
solution = Optim.optimize(obj_func,grad!,hess!,x0,method=NelderMead(),show_trace=false, iterations=10000)
println(solution)
x0 = Optim.minimizer(solution)
solution = Optim.optimize(obj_func,grad!,hess!,x0,method=Newton(),show_trace=false)
println(solution)
println(Optim.minimizer(solution))
obs_p = Optim.minimizer(solution)[1:N-1]
obs_t = Optim.minimizer(solution)[N:end]
#=
x0 = rand(2*(N-1))
@time solution = optimize(obj_func,grad!,hess!,x0,method=Newton(),show_trace=false)
println(solution)

x0 = rand(2*(N-1))
@time solution = optimize(obj_func,grad!,hess!,x0,method=Newton(),show_trace=false)
println(solution)
=#

# constrained problem
function obj_func(theta)
	if !all((0 .<= theta .<= 1)) # if not all elements in [0,1]
		return Inf
	else
		rho_vec = [1; theta[1:N-1]]
		t_vec = [t_lb; 0.0; theta[N:end] ; t_ub]
		profit = 0.0
		for i = 1:N-1
			f(t) = (rho_vec[i+1] - c)*d(rho_vec[i+1],t)
			int1 = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
			g(p) = d(p,t_vec[i+1])
			int2 = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
			inc = int1 + (1-t_cdf(t_vec[i+1]))*int2
			profit = profit + inc
		end
		return -profit
	end
end

function grad!(g_vec,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; 0.0; theta[N:end] ; t_ub]
	#price focs
	for i = 1:N-1 # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(rho_vec[i+1],t))
		res_rho = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(p,t_vec[i+1]))
		res_t = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		g_vec[i] = -res_rho
		g_vec[i+N-1-1] = -res_t
	end
end
function hess!(hess_mat,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; 0.0 ;theta[N:end] ; t_ub]
	#rho rho part of jac
	for i = 1:N-1 # i refers to scheduel part, not index position
		for j = 1:N-1
			if j == i
				f(t) = ((rho_vec[i+1] - c)*dd_dp_dp(rho_vec[i+1],t) + dd_dp(rho_vec[i+1],t) + (1-t_cdf(t))/(t_pdf(t))*(dd_dt_dp(rho_vec[i+1],t)))*t_pdf(t)
				res = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
				hess_mat[i,j] = -res
			else
				res = 0
				hess_mat[i,j] = -res
			end
		end
	end

	# rho t part of jac
	for i = 1:N-1
		for j = 2:N-1
			if j == i
				res = -((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1])/t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1],t_vec[i+1]))*t_pdf(t_vec[i+1])
				hess_mat[i,j+N-1-1] = -res
			elseif j == i+1
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1+1]) + (1 - t_cdf(t_vec[i+1+1]))/(t_pdf(t_vec[i+1+1]))*dd_dt(rho_vec[i+1],t_vec[i+1+1]))*t_pdf(t_vec[i+1+1])
				hess_mat[i,j+N-1-1] = -res
			else
				res = 0
				hess_mat[i,j+N-1-1] = -res
			end
		end
	end
	# t rho part of jac
	for i = 2:N-1
		for j = 1:N-1
			if j == i - 1
				res = ((rho_vec[i+1-1] - c)*dd_dp(rho_vec[i+1-1],t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]))*dd_dt(rho_vec[i+1-1],t_vec[i+1]))
				hess_mat[i+N-1-1,j] = -res
			elseif j == i
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(rho_vec[i+1],t_vec[i+1]))
				hess_mat[i+N-1-1,j] = res # no neg for some reason
			else
				res = 0
				hess_mat[i+N-1-1,j] = -res
			end
		end
	end
	# t t part of jac
	for i = 2:N-1
		for j = 2:N-1
			if j == i
				g(p) = (p - c)*dd_dp_dt(p,t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt_dt(p,t_vec[i+1]) + dd_dt(p,t_vec[i+1])*(-1.0)
				res = sparse_int(g,rho_vec[i+1], rho_vec[i+1-1])
				hess_mat[i+N-1-1,j+N-1-1] = -res
			else
				res = 0
				hess_mat[i+N-1-1,j+N-1-1] = -res
			end
		end
	end
end
#=
#checking grad/hess
testx = [3/5,1/5,1/5]
step = 1e-9
gtest1 = ones(3)
gtest2 = ones(3)
htest = ones(3,3)
eps = zeros(3)
step = 1e-9
eps[3]= step
upx = testx + eps
downx = testx - eps
est_grad = (obj_func(upx) - obj_func(downx))/(2*step)
println("Numerical grad: ", est_grad)
grad!(testx,gtest1)
println(gtest1)

grad!(upx,gtest1)
grad!(downx,gtest2)
hess!(testx,htest)
println("Numerical Hessian: ", (gtest1 - gtest2)/(2*step))
println(htest)
=#
x0 = rand(2*(N-1)-1)
@time solution = Optim.optimize(obj_func,grad!,hess!,x0,method=NewtonTrustRegion(),show_trace=false)
@time solution = Optim.optimize(obj_func,grad!,hess!,x0,method=NewtonTrustRegion(),show_trace=false)
@time solution = Optim.optimize(obj_func,grad!,hess!,x0,method=NewtonTrustRegion(),show_trace=false)
@time solution = Optim.optimize(obj_func,grad!,hess!,x0,method=NewtonTrustRegion(),show_trace=false)
println(solution)
println(Optim.minimizer(solution))
println("Profits: ",-obj_func(Optim.minimizer(solution)))

solution = Optim.optimize(obj_func,x0,method=NelderMead())
println(solution)
println(Optim.minimizer(solution))
println("Profits: ",-obj_func(Optim.minimizer(solution)))

solution = Optim.optimize(obj_func,x0,method=SimulatedAnnealing(),iterations=500000)
println(solution)
println(Optim.minimizer(solution))
println("Profits: ",-obj_func(Optim.minimizer(solution)))

## Identification test
#= the goal of this test is to take the observed price schedule and see if the parameters
can be recovered by solving the FOCs for 0 as a function of the params. =#
id_N = N
p(i,n) = 1 - i/(n - 1/2)
t(i,n) = (i - .5)/(n-.5)

p_vec = [1.0; [p(i,id_N) for i = 1:(id_N-1)]]
t_vec = [t_lb ; [t(i,id_N) for i = 1:(id_N-1)] ; t_ub]
println(p_vec)
println(t_vec)

function id_focs(theta)
	if (theta[2] < 1.0) | (theta[3] < 0.0) | (theta[1] < 0.0)
		res = Inf
	else
		ic = theta[1]
		i_pdf(x) = ks_dens(x,theta[2],theta[3])
		i_cdf(x) = ks_dist(x,theta[2],theta[3])
		rho_vec = [1.0; obs_p]
		t_vec = [t_lb; obs_t ; t_ub]
		#price focs
		res = 0.0
		for i = 1:id_N-1 # note index is shifted up by 1 to deal with 1 indexing
			f(t) = ((rho_vec[i+1] - ic)*dd_dp(rho_vec[i+1],t) + (1 - i_cdf(t))/i_pdf(t)*dd_dt(rho_vec[i+1],t))
			res_rho = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
			g(p) = ((p - ic)*dd_dp(p,t_vec[i+1]) + (1-i_cdf(t_vec[i+1]))/i_pdf(t_vec[i+1])*dd_dt(p,t_vec[i+1]))
			res_t = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
			res = res_rho^2 + res_t^2 + res
		end
	end
	return res
end

x0 = [0.0,1.5,1.5]
println(id_focs([c,b,a]))
solution = Optim.optimize(id_focs,x0,method=NelderMead(), g_tol=1e-12)
println(solution)
println(Optim.minimizer(solution))


#=
#### JuMP solution

#registering functions
function obj_func(params...)
	rho_vec = [1; collect(params[1:N-1])]
	t_vec = [t_lb; collect(params[N:end]);  t_ub]
	profit = 0.0
	for i = 1:N-1
		f(t) = (rho_vec[i+1] - c)*d(rho_vec[i+1],t)
		int1 = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = d(p,t_vec[i+1])
		int2 = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		inc = int1 + (1-t_cdf(t_vec[i+1]))*int2
		profit = profit + inc
	end
	return -profit
end

function grad!(g_vec,params...)
	rho_vec = [1; collect(params[1:N-1])]
	t_vec = [t_lb; collect(params[N:end]);  t_ub]

	#price focs
	for i = 1:N-1 # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t) + (1 - t_cdf(t))/t_pdf(t)*dd_dt(rho_vec[i+1],t))
		res_rho = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/t_pdf(t_vec[i+1])*dd_dt(p,t_vec[i+1]))
		res_t = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		g_vec[i] = -res_rho
		g_vec[i+N-1] = -res_t
	end
end

N = 3

#JuMP.register(:obj_func, 2,obj_func,grad!,hess!, autodiff=false)
JuMP.register(:obj_func, 2*(N-1),obj_func,grad!,autodiff=false)

m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
@variable(m,t1, start = 0.0)
@variable(m,t2, start = 1/3)
@variable(m,p1, start = 1/2)
@variable(m,p2, start = 1/3)
@NLconstraint(m,cons1,t1*(1-p1)==0)
@NLobjective(m,Min,obj_func(p1,p2,t1,t2))

# @time sol = solve(m)

println(getvalue(p1))
println(getvalue(p2))
println(getvalue(t1))
println(getvalue(t2))

#m2 = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
m2 = Model(solver=IpoptSolver(print_level = 0))
@variable(m2, s[1:(2*(N-1))])
x0 = [1/2,1/3,0.0,1/3]
for i = 1:2*(N-1)
	setvalue(s[i],x0[i])
end

@NLconstraint(m2,cons1,s[N]*(1-s[1])==0.0)
#setupperbound(s[N],0.0)
#setlowerbound(s[N],0.0)
#@NLobjective(m2,Min,obj_func(s[1],s[2],s[3],s[4]))
obj_str = "@NLobjective(m2,Min,obj_func("
for i = 1:2*(N-1)
	obj_str = obj_str*"s[$i],"
end
obj_str = obj_str[1:end-1]
obj_str = obj_str*"))"
eval(parse(obj_str))

@time sol2 = solve(m2)
println(getvalue(s))
#=
for i = 1:2*(N-1)
	setvalue(s[i],x0[i])
end
@time sol2 = solve(m2)
for i = 1:2*(N-1)
	setvalue(s[i],x0[i])
end
@time sol2 = solve(m2)
for i = 1:2*(N-1)
	setvalue(s[i],x0[i])
end
@time sol2 = solve(m2)
=#
=#
#=
#### MathProgBase
pert = 1e-12
function obj_func(theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	profit = 0.0
	for i = 1:N-1
		f(t) = (rho_vec[i+1] - c)*d(rho_vec[i+1],t)
		int1 = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = d(p,t_vec[i+1])
		int2 = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		inc = int1 + (1-t_cdf(t_vec[i+1]))*int2
		profit = profit + inc
	end
	return -profit
end

function grad!(g_vec,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	#price focs
	for i = 1:N-1 # note index is shifted up by 1 to deal with 1 indexing
		f(t) = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t) + (1 - t_cdf(t))/(t_pdf(t)+pert)*dd_dt(rho_vec[i+1],t))
		res_rho = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
		g(p) = ((p - c)*dd_dp(p,t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1])+pert)*dd_dt(p,t_vec[i+1]))
		res_t = sparse_int(g,rho_vec[i+1],rho_vec[i+1-1])
		g_vec[i] = -res_rho
		g_vec[i+N-1] = -res_t
	end
end
function hess!(hess_mat,theta)
	rho_vec = [1; theta[1:N-1]]
	t_vec = [t_lb; theta[N:end] ; t_ub]
	#rho rho part of jac
	for i = 1:N-1 # i refers to scheduel part, not index position
		for j = 1:N-1
			if j == i
				f(t) = ((rho_vec[i+1] - c)*dd_dp_dp(rho_vec[i+1],t) + dd_dp(rho_vec[i+1],t) + (1-t_cdf(t))/(t_pdf(t)+pert)*(dd_dt_dp(rho_vec[i+1],t)))*t_pdf(t)
				res = sparse_int(f,t_vec[i+1],t_vec[i+1+1])
				hess_mat[i,j] = -res
			else
				res = 0
				hess_mat[i,j] = -res
			end
		end
	end

	# rho t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				res = -((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1])/(t_pdf(t_vec[i+1]+pert)))*dd_dt(rho_vec[i+1],t_vec[i+1]))*t_pdf(t_vec[i+1])
				hess_mat[i,j+N-1] = -res
			elseif j == i+1
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1+1]) + (1 - t_cdf(t_vec[i+1+1]))/(t_pdf(t_vec[i+1+1])+pert)*dd_dt(rho_vec[i+1],t_vec[i+1+1]))*t_pdf(t_vec[i+1+1])
				hess_mat[i,j+N-1] = -res
			else
				res = 0
				hess_mat[i,j+N-1] = -res
			end
		end
	end
	# t rho part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i - 1
				res = ((rho_vec[i+1-1] - c)*dd_dp(rho_vec[i+1-1],t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]+pert))*dd_dt(rho_vec[i+1-1],t_vec[i+1]))
				hess_mat[i+N-1,j] = -res
			elseif j == i
				res = ((rho_vec[i+1] - c)*dd_dp(rho_vec[i+1],t_vec[i+1]) + (1-t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]+pert))*dd_dt(rho_vec[i+1],t_vec[i+1]))
				hess_mat[i+N-1,j] = res #not sure why no neg here
			else
				res = 0
				hess_mat[i+N-1,j] = -res
			end
		end
	end
	# t t part of jac
	for i = 1:N-1
		for j = 1:N-1
			if j == i
				g(p) = (p - c)*dd_dp_dt(p,t_vec[i+1]) + (1 - t_cdf(t_vec[i+1]))/(t_pdf(t_vec[i+1]+pert))*dd_dt_dt(p,t_vec[i+1]) + dd_dt(p,t_vec[i+1])*(-1.0)
				res = sparse_int(g,rho_vec[i+1], rho_vec[i+1-1])
				hess_mat[i+N-1,j+N-1] = -res
			else
				res = 0
				hess_mat[i+N-1,j+N-1] = -res
			end
		end
	end
end



#=
#checking grad/hess
testx = [3/5,1/5,1/5,3/5]
testx = [3/4,1/4,0.1,1/2]
step = 1e-9
gtest1 = ones(2*(N-1))
gtest2 = ones(2*(N-1))
htest = ones(2*(N-1),2*(N-1))
eps = zeros(2*(N-1))
step = 1e-9
eps[4]= step
upx = testx + eps
downx = testx - eps
est_grad = (obj_func(upx) - obj_func(downx))/(2*step)
println("Numerical grad: ", est_grad)
grad!(gtest1,testx)
println(gtest1)

grad!(gtest1,upx)
grad!(gtest2,downx)
hess!(htest,testx)
println("Numerical Hessian: ", (gtest1 - gtest2)/(2*step))
println(htest)
=#
#=
type WilsonNLP <: MathProgBase.AbstractNLPEvaluator
end

function MathProgBase.initialize(d::WilsonNLP, requested_features::Vector{Symbol})
	for feat in requested_features
		if !(feat in [:Grad, :Jac, :Hess])
			error("Unsupported feature $feat")
		end
	end
end

MathProgBase.features_available(d::WilsonNLP) = [:Grad, :Jac, :Hess]

MathProgBase.eval_f(d::WilsonNLP, x) = obj_func(x)

function MathProgBase.eval_g(d::WilsonNLP,g,x)
	g[1] = x[N]*(1-x[1])
end

function MathProgBase.eval_grad_f(d::WilsonNLP,grad_f,x)
	test_vec = similar(grad_f)
	grad!(test_vec,x)
	grad_f[:] = test_vec[:]
end

MathProgBase.jac_structure(d::WilsonNLP) = [1,1,1,1], [1,2,3,4] # should vary with N
MathProgBase.hesslag_structure(d::WilsonNLP) = [1,2,2,3,3,3,4,4,4,4],[1,1,2,1,2,3,1,2,3,4] # should vary with N

function MathProgBase.eval_jac_g(d::WilsonNLP,J,x)
	J[1] = -x[N]
	J[2] = 0.0
	J[3] = (1-x[1])
	J[4] = 0.0
end

function MathProgBase.eval_hesslag(d::WilsonNLP, H, x, σ, μ)
	# Only specifying lower left triangle
	test_mat = ones(4,4)
	hess!(test_mat,x)
	H[1] = σ * test_mat[1,1] #1,1
	H[2] = σ * test_mat[2,1] #2,1
	H[3] = σ * test_mat[2,2] #2,2
	H[4] = σ * test_mat[3,1] #3,1
	H[5] = σ * test_mat[3,2] #3,2
	H[6] = σ * test_mat[3,3] #3.3
	H[7] = σ * test_mat[4,1] #4,1
	H[8] = σ * test_mat[4,2] #4,2
	H[9] = σ * test_mat[4,3] #4,3
	H[10] = σ * test_mat[4,4] #4,4

	# Constraint
	H[4] += μ[1]*(-1) # lower triangle only
end
function num_sol(solver=IpoptSolver(print_level=0))
	m = MathProgBase.NonlinearModel(solver)
	l = zeros(4)
	u = ones(4)
	lb = [0.0]
	ub = [0.0]

	MathProgBase.loadproblem!(m,4,1,l,u,lb,ub,:Min,WilsonNLP())
	start = [1/2,1/3,0.1,1/3]
	MathProgBase.setwarmstart!(m,start)

	MathProgBase.optimize!(m)
	x = MathProgBase.getsolution(m)
	println(x)
end

num_sol()

=# =#
