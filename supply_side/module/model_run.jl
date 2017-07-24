include("LiquorQD.jl")
using LiquorQD, DataFrames, DataArrays

#= Set up of coefficients. Product and group utility goes in product b/c those
don't multiply out well. The remainder go in a DemandCoefs object with an associated
weight =#
const group_utils = [-.9980637;.8108195;-.2365993;-.7477314] # order is gin,vod,rum,sch
const prod_utils= [-3.121894;-3.121704;-4.801027;-5.014518;-4.846033;-3.435482;
-4.777261;-2.949315;-4.601557;-4.360644;-5.3786;-5.390806;-3.930106;-4.342176;
-4.342773;-4.797522;-5.466756;-4.916693;-3.724147;-5.618995;-4.223;-3.830568;
-4.448236;-4.904963;-6.201783;-3.911629;-4.470888;-6.215972;-5.339201;-4.433005;
-5.04312;-5.455826;-4.350398;-4.999972;-5.11608;-5.068243;-5.927884;-5.143831;
-4.432692;-5.631095;-4.801862;-5.574736;-5.192336;-4.960927;-4.160068;-3.903724;
-5.181415;-3.886327;-5.344421;-5.285977;-4.1023;-5.637782;-4.822594;-3.958131;
-6.565579;-4.548577;-5.601947;-3.987991;-4.912977;-4.840165;-6.880695;-5.605676;
-5.188475;-5.377887;-4.09764;-5.108016;-5.710825;-5.113058;-5.991073;-5.567423;
-6.625885;-5.299186;-6.272768;-5.579566;-5.928923;-4.585755;-5.869962;-5.508581;
-4.829222;-6.51405;-5.81735;-6.198225;-6.39506;-5.299936;-7.157789;-5.640015;
-5.26857;-4.868245;-4.543375;-4.700213;-4.514191;-5.530586;-4.470521;-4.471254;
-5.435326;-4.380613;-6.817262;-6.6457;-6.61709;-6.495391;-.1659607]

coefs_lt30 = DemandCoefs(-.1288769,-.1549591,.0300896)
coefs_30_to_70 = DemandCoefs(-.0373797,-.3626817,.0154656)
coefs_70_to_100 = DemandCoefs(-.0312302,-.2814536,.0191501)
coefs_100_plus = DemandCoefs(-.028152,-.2519711, -.0019593) # price coef is estimated positive. causes problems

const inc_weights = [.2055;.3970;.2059;.1917]

const coef_array = [coefs_lt30;coefs_30_to_70;coefs_70_to_100;coefs_100_plus]


#= Initial set up of product data. Need to:
  1) import data
  2) initilize a product object for each product in data
  3) Construct markets out of products
=#

prod_char_df = readtable("E:/Box/Box Sync/Snuffles Backup/gh8728/Projects/Liquor/Data/combined/merged_sample.csv")
# estimated price schedule
rho_cols = [:actual_p0, :actual_p1, :actual_p2, :actual_p3, :actual_p4, :actual_p5, :actual_p6,
  :actual_p7, :actual_p8, :actual_p9]
cutoff_q_cols = [:disc_q0, :disc_q1, :disc_q2, :disc_q3, :disc_q4, :disc_q5, :disc_q6,
  :disc_q7, :disc_q8, :disc_q9]
mkt_ids = Array(levels(prod_char_df[:,:mkt])) # get market ids
prod_array = Liquor[]
markets_array = Market[]
for m in mkt_ids
  mkt_df = prod_char_df[(prod_char_df[:,:mkt] .== m),:] # get products in market
  mkt_vec = Liquor[]
  mkt_ps_lookup = Dict{Liquor,Tuple{Array{Float64,1},Array{Float64,1}}}() # dictionary with keys as products and tuple of prices and ffs as values.
  for row in eachrow(mkt_df)
    # constructing product objects
    id = row[:product]
    group_d_vec = [row[:d_g_gin],row[:d_g_vod],row[:d_g_rum],row[:d_g_sch]] # get producrt group dummy values
    g_util = sum(group_d_vec.*group_utils) # get group utility
    prod = Liquor(id,row[:price],row[:price],row[:imported],row[:proof],g_util,prod_utils[id],nothing)
    push!(mkt_vec,prod) # push to market
    push!(prod_array,prod) # push to total product array
    # extracting data needed to construct price schedule objects
    matched_ps = (row[:_merge_ps] == "matched (3)")
    case = (row[:discount_size_type] == "CASE")
    if matched_ps
      obs_rhos = Float64[]
      for col in rho_cols
        if !isna(row[col])
          push!(obs_rhos,row[col])
        end
      end
      obs_cutoff_q = Float64[]
      for col in cutoff_q_cols
        if !isna(row[col])
          push!(obs_cutoff_q,row[col])
        end
      end
      if case
        obs_rhos = obs_rhos ./ 12.0
        obs_cutoff_q = obs_cutoff_q .* 10.0
      end
      obs_ff = [0.0] # first fixed fee is by definition, and corresponds to A1
      for k = 2:length(obs_cutoff_q)
  			diff = (obs_cutoff_q[k] - 1)*(obs_rhos[k-1] - obs_rhos[k])
  			res = diff + obs_ff[k-1]
  			push!(obs_ff,res)
  		end
      mkt_ps_lookup[prod]=(obs_rhos,obs_ff)
    end
  end
  date_str = mkt_df[1,:date_m]
  mkt_y = parse(Int64,date_str[1:4])
  mkt_m = parse(Int64,date_str[6:end])
  tmp_mkt = Market(mkt_vec,mkt_m,mkt_y)
  # constructing price schedule objects
    for j in tmp_mkt.products
      try
        tmp_rhos = mkt_ps_lookup[j][1]
        tmp_ff = mkt_ps_lookup[j][2]
        tmp_t = obs_lambdas(tmp_rhos,tmp_ff,j,coef_array,inc_weights,tmp_mkt,10000000.0) # M probably not right. Need to pin down a better number
        tmp_ps = PriceSched(tmp_rhos,tmp_t,(length(tmp_rhos)+1))
        j.ps = tmp_ps # Setting price schedule for product
      end
    end
  push!(markets_array,tmp_mkt)
end

for m in markets_array[1:1]
  for j in m.products
    if !isnull(j.ps)
      println("working with product", j)
      tmp_ps = get(j.ps) # becase the ps field is nullable, need to use get
      dev_ps = dev_gen(tmp_ps,0.05)
      dev_ps = dev_ps[rand(1:end,200)] # 200 random ineqaulities. Speeds up computation
      print("Pre-calculating retail prices. ")
      pre_calc = Dict{Int64,Float64}[]
      for s in dev_ps
        tmp_dict = Dict(i => p_star(s.rhos[i],j,coef_array,inc_weights,m) for i in 1:(s.N-1))
        push!(pre_calc,tmp_dict)
      end
      println("Done.")
      min_rho = minimum(tmp_ps.rhos)
      sol,xtrace,ftrace = optimize_moment(tmp_ps,dev_ps,j,coef_array,inc_weights,m,25000,pre_calc,x0=[min_rho,1.0,1.0])
      println(sol)
      trace = [vcat(xtrace'...) ftrace]
      out_str = "traces/trace_"*string(m.year)*"_"*string(m.month)*"_"*string(j.id)*".csv"
      writedlm(out_str,trace)
    end
  end
end

#=
m = markets_array[1]
j = m.products[1]
tmp_ps = get(j.ps)
println("Generating deviation price schedules.")
dev_ps = dev_gen(tmp_ps,0.2)
println("Pre-calculating retail prices.")
pre_calc = Dict{Int64,Float64}[]
for s in dev_ps
  tmp_dict = Dict(i => p_star(s.rhos[i],j,coef_array,inc_weights,m) for i in 1:(s.N-1))
  push!(pre_calc,tmp_dict)
end

println("b, Q")
for b = 0.5:.5:20
  tmp_params = WholesaleParams(8.0,1.0,b)
  res = moment_obj_func(tmp_ps,dev_ps,tmp_params,j,coef_array,inc_weights,m,pre_calc)
  println(b,",",res,"\n")
end
=#
