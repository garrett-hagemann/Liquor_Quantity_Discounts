# addprocs(68)
# addprocs(2)

data_file_arg = ARGS[1] # needs a command line argument passed

@everywhere srand(69510606) #seeding random number gen

@everywhere include("LiquorQD.jl")
@everywhere include("DataSetup.jl")
@everywhere using DataSetup, LiquorQD, DataFrames, DataArrays

# setting up data
markets_array = data_setup(data_file_arg)

# importing estimated coefficients
ss_coefs_df = readtable("trace_sum.csv")
@everywhere est_coefs = Dict{Tuple{Int64,Int64,Int64},Tuple{Float64,Float64}}() # epty dict. Keys are (y,m,p) and value is (c,b)
for row in eachrow(ss_coefs_df)
    key = (row[:year],row[:month],row[:prod])
    val = (row[:c_mode],row[:b_mode])
    est_coefs[key]=val
end

@everywhere function cf_calcs(tup::Tuple{Market,Liquor},coefs_dict::Dict{Tuple{Int64,Int64,Int64},Tuple{Float64,Float64}})
    max_seg = 15
    m = tup[1]
    j = tup[2]
    out = Dict{Tuple{Market,Liquor},Array{Float64,1}}()
    key = (m.year,m.month,j.id)
    println("working with product $(j.id) in $(m.year)-$(m.month)")
    if key in keys(coefs_dict)
        cf_M = 2000000.0
        wp = WholesaleParams(coefs_dict[key][1],1.0,coefs_dict[key][2])
        try
            obs_ps = get(j.ps) # because it's nullable
            println("Observed N: ",obs_ps.N)
            println("Estimated parameters: ", wp)

            ps_dict = optimal_price_sched(wp,max_seg,j,nested_coefs,obs_inc_dist,m) # calculating all necessary price schedules

            # convergence to max profit check
            prof_diff_rel = (wholesaler_profit(ps_dict[max_seg],wp,j,nested_coefs,obs_inc_dist,m) - wholesaler_profit(ps_dict[max_seg-1],wp,j,nested_coefs,obs_inc_dist,m))/wholesaler_profit(ps_dict[max_seg-1],wp,j,nested_coefs,obs_inc_dist,m)
            println("Relative difference in profit between $max_seg and $(max_seg-1) segments: ", prof_diff_rel)

            base_ps = ps_dict[obs_ps.N] # this ensures inequality defining zeta is true. Can't use observed
            ps_up = ps_dict[obs_ps.N+1]
            ps_dn = ps_dict[obs_ps.N-1]
            ps_max = ps_dict[max_seg]

            #finding zeta
            base_w_profit =  wholesaler_profit(base_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            profit_up = wholesaler_profit(ps_up,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            profit_dn = wholesaler_profit(ps_dn,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            profit_max = wholesaler_profit(ps_max,wp,j,nested_coefs,obs_inc_dist,m)*cf_M

            print("Recovering zeta. ")
            zeta_lb = profit_up - base_w_profit
            zeta_ub = base_w_profit - profit_dn
            zeta_mid = (zeta_lb + zeta_ub)/2.0
            println("Bounds on zeta: [",zeta_lb,",",zeta_ub,"].")

            # change in  wholesaler profit
            lin_ps = ps_dict[2]
            lin_w_profit = wholesaler_profit(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            delta_w_profit = (lin_w_profit - 2*zeta_mid) - (base_w_profit - zeta_mid*obs_ps.N)
            println("Change in wholesaler profits under linear price: ", delta_w_profit)

            # change in retailer profit
            base_r_profit = ps_ret_profits(base_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            lin_r_profit = ps_ret_profits(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            delta_r_profit = lin_r_profit - base_r_profit
            println("Change in retailer profit: ", delta_r_profit)

            # change in consumer surplus
            base_cs = ps_cons_surplus(base_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            lin_cs = ps_cons_surplus(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)*cf_M
            delta_cs = lin_cs - base_cs
            println("Change in consumer surplus: ", delta_cs)

            base_avg_p = ps_avg_ret_price(base_ps,wp,j,nested_coefs,obs_inc_dist,m)
            lin_avg_p = ps_avg_ret_price(lin_ps,wp,j,nested_coefs,obs_inc_dist,m)
            println("Avg retail price under observed schedule: ", base_avg_p)
            println("Avg retail price under linear price: ", lin_avg_p)

            base_wavg_p = ps_swavg_ret_price(base_ps,wp,j,nested_coefs,obs_inc_dist,m)
            println("Weighted avg retail price under observed schedule: ", base_wavg_p)

            #= calculating required points to find difference in profit as a function of
            retailer type =#
            t_r_prof_array = Float64[]
            lin_r_p_star = p_star(lin_ps.rhos[1],j,nested_coefs,obs_inc_dist,m)
            lin_r_share = share(lin_r_p_star,j,nested_coefs,obs_inc_dist,m)
            base_ff = recover_ff(base_ps,j,nested_coefs,obs_inc_dist,m)

            for (rho,t,ff) in zip(base_ps.rhos,base_ps.t_cuts,base_ff)
                tmp_p_star = p_star(rho,j,nested_coefs,obs_inc_dist,m)
                tmp_t_r_prof = t*cf_M*((lin_r_p_star - lin_ps.rhos[1])*lin_r_share - (tmp_p_star - rho)*share(tmp_p_star,j,nested_coefs,obs_inc_dist,m)) + ff*cf_M
                push!(t_r_prof_array,tmp_t_r_prof)
            end
            # adding one more for type = 1
	    last_t_p_star = p_star(base_ps.rhos[end],j,nested_coefs,obs_inc_dist,m)
            last_t_r_prof = 1.0*cf_M*((lin_r_p_star - lin_ps.rhos[1])*lin_r_share - (last_t_p_star - base_ps.rhos[end])*share(last_t_p_star,j,nested_coefs,obs_inc_dist,m)) - base_ff[end]*cf_M
            push!(t_r_prof_array,last_t_r_prof)

            res_array = [zeta_lb; zeta_ub; zeta_mid; (base_w_profit) ; (lin_w_profit); profit_max; 
                base_r_profit ; lin_r_profit ; base_cs ; lin_cs ; base_avg_p ; lin_avg_p ; base_wavg_p ; base_ps.t_cuts; t_r_prof_array]
            out[(m,j)] = res_array
        end
    else
        println("Product has no parameter estimates.")
    end
    return out
end

mkts_for_est = [(m,j) for m in markets_array for j in m.products] # array of tuples
res = pmap((x)->cf_calcs(x,est_coefs),mkts_for_est)

#outfile = open("counter_fact_calcs.csv","w")
#write(outfile,"year,month,prod,zeta_lb,zeta_ub,zeta_mid,delta_w_profit,delta_r_profit,delta_cs,base_avg_p,lin_avg_p,base_wavg_p","\n")
col_heads = "year,month,prod,zeta_lb,zeta_ub,zeta_mid,base_w_profit,lin_w_profit,max_w_profit,base_r_profit,lin_r_profit,base_cs,lin_cs,base_avg_p,lin_avg_p,base_wavg_p"
for d in res
    for (key,v) in d # should only be single element, but just in case
        m = key[1]
        j = key[2]
	tmp_ps = get(j.ps)
        type_string = join(["t_cut_$i" for i=1:(tmp_ps.N-1)],",")
        r_prof_string = join(["r_t_prof_$i" for i=1:(tmp_ps.N-1)],",")

        row = [m.year;m.month;j.id;v]'
        out_str = "cf/cf_$(m.year)_$(m.month)_$(j.id).csv"
	header = join([col_heads;type_string;r_prof_string],",")
        f_out = open(out_str,"w")
        write(f_out,header,"\n")
        writecsv(f_out,row)
        close(f_out)
    end
end
