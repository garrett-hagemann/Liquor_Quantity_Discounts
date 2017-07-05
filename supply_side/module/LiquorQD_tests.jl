include("LiquorQD.jl")
using LiquorQD

# setting up fake products, markets, and coefficients
test_prod1 = Liquor(1,13.56,13.56,0,90.0,0.0,0.0)
test_prod2 = Liquor(2,9.30,9.30,0,40.0,0.0,0.0)

test_mkt = Market([test_prod1],6,2009)

test_coefs = DemandCoefs(-0.4,0.0,0.0)
test_weights = [1.0]

# testing density functions and integral
ks_a = 2.0
ks_b = 2.0

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
test_w_params = WholesaleParams(0.0,1.0,1.0,4)
test_ps = optimal_price_sched(test_w_params,test_prod1,[test_coefs],test_weights,test_mkt)
println("Optimal price schedule: ", test_ps)
println("Profit at optimal schedule: ", wholesaler_profit(test_ps,test_w_params,test_prod1,[test_coefs],test_weights,test_mkt))
