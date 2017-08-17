# addprocs(68)
#addprocs(2)

@everywhere srand(69510606) #seeding random number gen

@everywhere include("LiquorQD.jl")
@everywhere include("DataSetup.jl")
@everywhere using DataSetup, LiquorQD, DataFrames, DataArrays

# setting up data
markets_array = data_setup()

# importing estimated coefficients
ss_coefs_df = readtable("traces/trace_sum.csv")
coefs_dict = Dict{Tuple{Int64,Int64,Int64},Tuple{Float64,Float64}}() # epty dict. Keys are (y,m,p) and value is (c,b)
for row in eachrow(ss_coefs_df)
    key = (row[:year],row[:month],row[:prod])
    val = (row[:c_mode],row[:b_mode])
    coefs_dict[key]=val
end

for m in markets_array
    for j in m.products
        key = (m.year,m.month,j.id)
        if key in keys(coefs_dict)
            print("Recovering zeta for product $(j.id) in $(m.year)-$(m.month). ")
            wp = WholesaleParams(coefs_dict[key][1],1.0,coefs_dict[key][2])
            tmp_ps = get(j.ps) # because it's nullable
            #finding zeta
            @time base_ps = optimal_price_sched(wp,(tmp_ps.N),j,nested_coefs,obs_inc_dist,m) # this ensures inequality defining zeta is true. Can't use observed
            base_w_profit =  wholesaler_profit(tmp_ps,wp,j,nested_coefs,obs_inc_dist,m)
            @time ps_up = optimal_price_sched(wp,(tmp_ps.N+1),j,nested_coefs,obs_inc_dist,m)
            profit_up = wholesaler_profit(ps_up,wp,j,nested_coefs,obs_inc_dist,m)
            @time ps_dn = optimal_price_sched(wp,(tmp_ps.N-1),j,nested_coefs,obs_inc_dist,m)
            profit_dn = wholesaler_profit(ps_dn,wp,j,nested_coefs,obs_inc_dist,m)

            println(base_ps)
            println(ps_up)
            println(ps_dn)
            println(base_w_profit)
            println(profit_up)
            println(profit_dn)


            zeta_lb = profit_up - base_w_profit
            zeta_ub = base_w_profit - profit_dn
            zeta_mid = (zeta_lb + zeta_ub)/2.0
            println("Bounds on zeta: [",zeta_lb,",",zeta_ub,"].")

            # change in  wholesaler profit
            lin_ps = optimal_price_sched(wp,2,j,nested_coefs,obs_inc_dist,m)
            lin_w_profit = wholesaler_profit(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            delta_w_profit = (lin_w_profit - 2*zeta_mid) - (base_w_profit - zeta_mid*tmp_ps.N)
            println("Change in wholesaler profits under linear price: ", delta_w_profit)

            # change in retailer profit
            base_r_profit = ps_ret_profits(tmp_ps,wp,j,nested_coefs,obs_inc_dist,m)
            lin_r_profit = ps_ret_profits(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            delta_r_profit = lin_r_profit - base_r_profit
            println("Change in retailer profit: ", delta_r_profit)

            # change in consumer surplus
            base_cs = ps_cons_surplus(tmp_ps,wp,j,nested_coefs,obs_inc_dist,m)
            lin_cs = ps_cons_surplus(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            delta_cs = lin_cs - base_cs
            println("Change in consumer surplus: ", delta_cs)

            base_avg_p = ps_avg_ret_price(tmp_ps,wp,j,nested_coefs,obs_inc_dist,m)
            lin_avg_p = ps_avg_ret_price(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            println("Avg retail price under observed schedule: ", test_avg_p)
            println("Avg retail price under linear price: ", lin_avg_p)

            base_wavg_p = ps_swavg_ret_price(tmp_ps,wp,j,nested_coefs,obs_inc_dist,m)
            lin_wavg_p = ps_swavg_ret_price(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            println("Weighted avg retail price under observed schedule: ", test_wavg_p)
            println("Weighted avg retail price under linear price: ", lin_wavg_p)


        end
    end
end
