include("LiquorQD.jl")
using LiquorQD, DataFrames, DataArrays

#= Set up of coefficients. Product and group utility goes in product b/c those
don't multiply out well. The remainder go in a DemandCoefs object with an associated
weight =#
const size_utils = [.151699,.6506976,-.1413254] # order is 750ml, 1L, 1.75L
const prod_utils= [.9984359; 2.222908; -.2574402; 1.29754; .8887798; -.9459865; -.9818617;
-.974148; 1.211948; -.6472216; 1.93983; 1.943144; -.5035014; -.7443478;
-.8400436; -.5358744; 1.673104; 1.615472; -1.113532; -1.441404; 1.677565;
-.0269081; -1.235775; -.0093349; 2.235452; 1.241808; -.103252; -.0133208;
-.537227; 1.561229; 1.520777; .6162316; 1.571745; .1175427; -.776834; 2.242108;
-.3264315; -1.750993; -.4174726; .9231605; 1.797444; -1.18543; 1.98275; -.0485804;
-2.308151; .6849962; -2.050701; .8454664; 1.716071; -1.960274; -.4715936; .0917088;
1.788842; -1.454759; -1.672869; 1.066978; .0430603; 1.006194; -2.148041; -1.867818;
.7396636; .5350933; -1.569208; 1.195685; .3580877; 1.537322; -1.00704; 1.750797;
1.117172; -1.534264; 1.363423; 1.82211; -1.225709; .9798436; 2.132275; 1.892909;
1.866594; -.6740343; -2.048624; -2.685772; 1.213439; -1.622729; 1.026525; -1.385786;
-1.592726; -.3203502; 1.124569; .6860642; .4910821; .0602244; .6011467; .989683;
-2.809397; 1.772477; -2.359167; -1.26707; -.722626; -1.295764; 1.800456; 1.337842;
.3584706; 1.100851; .3543156; -.2078676; .7758802; -.2940698; 2.612169; 1.86315]

const nested_coefs = DemandCoefs(-.224743, # price
                            -.8281927, # imported
                             -.0075732, # proof
                             0.0, #.0202472, # price*y // setting to zero because it causes positive respones
                             .0780103, # imported*y
                              .0002535, # proof*y
                              Dict("oo" => 1,
                              "gin" => 1.255592,
                              "vod" => 1.38928,
                              "rum" => .555316,
                              "sch" => 1.264724,
                              "brb" => 1.822675,
                              "whs" => .2343492,
                              "teq" => 1.671329,
                              "otr" => .3359661))

inc_levels = Float64[2500; 6500; 9000; 11000; 13000; 17500; 22500; 27500; 32500; 37500; 42500; 47500; 55000; 65000; 85000; 150000]
inc_weights = [0.34; 2.36; 1.21; 2.63; 3.91; 2.02; 4.99; 9.84; 2.76; 2.83; 4.45; 2.49; 11.73; 4.38; 24.14; 19.89]

const obs_inc_dist = IncomeDist(inc_levels,inc_weights)
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
    size_d_vec = [row[:d_s_750ML],row[:d_s_1L],row[:d_s_175L]] # get size dummy values
    s_util = sum(size_d_vec.*size_utils) # get group utility
    prod = Liquor(id,row[:_type],row[:avg_price_inter],row[:avg_price_inter],row[:imported],row[:proof],s_util,prod_utils[id],nothing)
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
        tmp_t = obs_lambdas(tmp_rhos,tmp_ff,j,nested_coefs,obs_inc_dist,tmp_mkt,10000000.0) # M probably not right. Need to pin down a better number
        tmp_ps = PriceSched(tmp_rhos,tmp_t,(length(tmp_rhos)+1))
        j.ps = tmp_ps # Setting price schedule for product
      end
    end
  push!(markets_array,tmp_mkt)
end

for m in markets_array[1:1]
  for j in m.products
    println("test: ", j)
    if !isnull(j.ps)
      println("working with product", j)
      tmp_ps = get(j.ps) # becase the ps field is nullable, need to use get
      dev_ps = dev_gen(tmp_ps,0.05)
      dev_ps = dev_ps[rand(1:end,200)] # 200 random ineqaulities. Speeds up computation
      print("Pre-calculating retail prices. ")
      pre_calc = Dict{Int64,Float64}[]
      for s in dev_ps
        tmp_dict = Dict(i => p_star(s.rhos[i],j,nested_coefs,obs_inc_dist,m) for i in 1:(s.N-1))
        push!(pre_calc,tmp_dict)
      end
      println("Done.")
      min_rho = minimum(tmp_ps.rhos)
      sol,xtrace,ftrace = optimize_moment(tmp_ps,dev_ps,j,nested_coefs,obs_inc_dist,m,25000,pre_calc,x0=[min_rho,1.0,1.0])
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
