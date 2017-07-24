include("LiquorQD.jl")
using LiquorQD
srand(69510606) #seeding random number gen

# setting up fake products, markets, and coefficients
test_prod1 = Liquor(1,13.56,13.56,0,90.0,0.0,0.0,nothing) # null price sched
test_prod2 = Liquor(2,9.30,9.30,0,40.0,0.0,0.0,nothing) # null price sched

test_mkt = Market([test_prod1],6,2009)

test_coefs = DemandCoefs(-0.12,0.0,0.0)
test_weights = [1.0]

# testing density functions and integral
ks_a = 1.0
ks_b = 1.0

println("Integrating density: ", sparse_int(x->ks_dens(x,ks_a,ks_b),0.0,1.0)) # should be 1
println("F(.5): ", sparse_int(x->ks_dens(x,ks_a,ks_b),0.0,0.5)) # should be 1/2 when a=b=1.0
println("Closed form F(.5): ", ks_dist(0.5,ks_a,ks_b))

# testing share function
test_price = 13.56
res = share(test_price,test_prod1,[test_coefs],test_weights,test_mkt)
println("Share at $test_price: ", res)

#testing p_star function
test_mc = 10.0
ps_res = p_star(test_mc,test_prod1,[test_coefs],test_weights,test_mkt)
println("Optimal price at $test_mc: ", ps_res)

# Testing price schedule optimization
test_w_params = WholesaleParams(2.0,2.0,2.0)
test_N = 6
test_ps = optimal_price_sched(test_w_params,test_N,test_prod1,[test_coefs],test_weights,test_mkt)
println("Optimal price schedule: ", test_ps)
println("Profit at optimal schedule: ", wholesaler_profit(test_ps,test_w_params,test_prod1,[test_coefs],test_weights,test_mkt))
# testing deviation generation
test_δ = 0.025
test_devs = dev_gen(test_ps,test_δ)
test_devs = test_devs[rand(1:end,200)] # N random ineqaulities

# pre-calculating optimal retail prices since they don't change with wholesaler params
test_pre_calc = Dict{Int64,Float64}[]
for s in test_devs
  tmp_dict = Dict(i => p_star(s.rhos[i],test_prod1,[test_coefs],test_weights,test_mkt) for i in 1:(s.N-1))
  push!(test_pre_calc,tmp_dict)
end

# testing moment obj function evaluation
test_w_params2 = WholesaleParams(0.0,3.0,1.5)
test_Q = moment_obj_func(test_ps,test_devs,test_w_params,test_prod1,[test_coefs],test_weights,test_mkt,test_pre_calc)
println("Q(true params) = ", test_Q) # should b 0
test_Q2 = moment_obj_func(test_ps,test_devs,test_w_params2,test_prod1,[test_coefs],test_weights,test_mkt,test_pre_calc)
println("Q(other params) = ", test_Q2) # should be > 0

# testing moment optimization to recover params
optimize_moment(test_ps,test_devs,test_prod1,[test_coefs],test_weights,test_mkt,1,test_pre_calc,x0=[0.0,1.0,1.0])
@time test_recovered_params,test_xtrace,test_ftrace = optimize_moment(test_ps,test_devs,test_prod1,[test_coefs],test_weights,test_mkt,50000,test_pre_calc,x0=[0.0,1.0,1.0])
println("Recovered Parameters: ", test_recovered_params)
test_trace = [vcat(test_xtrace'...) test_ftrace]
writedlm("test_dens.csv",test_trace)

# mapping shape of objective function
#=
println("b, Q")
for b = 0.5:.5:20
  tmp_params = WholesaleParams(b,6.0,6.0)
  res = moment_obj_func(test_ps,test_devs,tmp_params,test_prod1,[test_coefs],test_weights,test_mkt,test_pre_calc)
  print(b,",",res,"\n")
end
=#
