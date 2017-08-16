module DataSetup

include("LiquorQD.jl")
using LiquorQD, DataFrames, DataArrays

const size_utils = [0.0,0.0,0.0] # order is 750ml, 1L, 1.75L
const prod_utils= [-1.031654;  -0.5986714;  2.719251;  -2.136661;  -0.7166995;  -2.262;  -2.120311;  -1.541387;
-2.646299;  -2.654972;  -1.434078;  -0.0966979;  -2.604326;  -2.413193;  -2.50415;  -2.278296;
-2.025461;  -2.939136;  -2.366031;  -2.627487;  -2.73191;  -0.7602613;  -1.76746;  -0.9681924;
-2.647136;  -3.143549;  -3.056831;  -1.480597;  -2.202042;  -0.7426222;  -2.521247;  -0.3321612;
-1.883391;  -2.37551;  -0.0969018;  -2.71941;  -0.0398781;  -2.479553;  -3.223033;  -4.171165;
-1.537278;  -3.579824;  -2.600178;  -3.214601;  -2.201891;  -2.242283;  -2.909715;  -3.208443;
-3.047639;  -4.36851;  -4.14869;  -4.302991;  -3.446561;  -0.4603544;  -1.889882;  -2.320702;
-3.668417;  -1.162365;  -3.543497;  -3.725587;  -4.398088;  -2.799301;  -2.832664;  -3.511715;
-2.207737;  -2.264005;  -2.970694;  -2.267005;  -5.102547;  -4.635895;  -3.291156;  -4.526457;
-2.711992;  -1.847489;  -3.467429;  -3.002993;  -2.899891;  2.297307;  -3.361088;  -2.685769;
-2.579266;  -1.914592;  -4.353611;  -3.701298;  -3.726989;  -3.745329;  -3.262226;  -0.8323384;
-1.837014;  -3.906951;  -0.7877522;  -1.406462;  -1.866976;  -1.028922;  -1.593316;  -3.226181;
-0.8440271;  -4.112507;  -3.948755;  -4.13223;  -3.119976;  -2.607345;  -3.573856;  -3.489196;
-4.073018;  -0.4705465;  -2.96645;  -2.777133;  -4.992293;  -3.35968;  -2.760934;  -4.169735;
-5.476136;  -3.135815;  -3.467021;  -2.619442;  -0.7818482;  -4.255969;  -4.730402;  -3.705244;
-4.757192;  -2.710783;  -4.723884;  -3.884752;  -0.5518308;  -2.169723;  -3.857576;  -4.087936;
-5.246792;  -1.449692;  -4.304432;  -4.742753;  -3.562025;  -5.305832;  -5.286488;  -5.105905;
-3.918621;  -3.978658;  -3.249214;  -4.694763;  -3.242497;  -4.142702;  -3.063908;  -4.552217;
-0.8854282;  -3.007056;  -4.586837;  -4.343445;  -1.675508;  -5.384725;  -4.351588;  -4.493908;
-1.829689;  -3.449844;  -3.834611;  -3.630683;  -4.76543;  -4.508183;  -4.113295;  2.337819;
-5.7323;  -4.490103;  -4.951924;  -0.9476656;  -4.890666;  -3.23345;  -3.968342;  -5.624683;
-2.367878;  -4.28221;  -4.702244;  -1.621754;  -2.360367;  -3.128938;  -5.312097;  -5.652628;
-5.148921;  -4.865083;  -4.679488;  -4.805693;  -5.244291;  -3.566109;  -3.501546;  -5.579437;
-2.110095;  -3.881827;  -3.088027;  -5.167253;  -3.314906;  -3.843971;  -1.675156;  -2.342578;
-3.314817;  -0.9385336;  -4.097188;  -1.996012;  -5.252544;  -4.562477;  -3.46129;  -4.061216;
-5.464092;-4.311072;-3.505353;-4.919189;-3.002201;-4.591751;-3.96926;-5.012443;-3.355068;-5.073336;
-2.593191;-4.422889;-4.37292;-0.7637211;-3.551962;-3.01347;-4.800927;-3.354938;-5.143094;-4.038786;
-3.781936;-5.678621;-6.073475;-4.826704;-4.586513;-4.349475;-6.103028;-4.223599;-4.233743;-1.581903;
-4.372692;-3.39948;-2.542516;-3.287263;-6.033693;-4.714857;-5.936033;-5.880964;-4.488038;-4.623158;
-5.179918;-4.234136;-4.884122; -1.525883; -5.316571; -3.39473; -4.659103; -5.109047; -3.767755; -2.927771;
0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0]

const nested_coefs = DemandCoefs(-.1596846, # price
                            0.0, # imported // unused currently
                             0.0, # proof // unused currently
                             .0079275, # price*y
                             0.0, # imported*y // unused currently
                              0.0, # proof*y // unused currently
                              Dict("oo" => 1.0, # nesting params // unused currently
                              "gin" => 1.0,
                              "vod" => 1.0,
                              "rum" => 1.0,
                              "sch" => 1.0,
                              "brb" => 1.0,
                              "whs" => 1.0,
                              "teq" => 1.0,
                              "otr" => 1.0))

inc_levels = Float64[2500; 6500; 9000; 11000; 13000; 17500; 22500; 27500; 32500; 37500; 42500; 47500; 55000; 65000; 85000; 150000]
inc_weights = [0.34; 2.36; 1.21; 2.63; 3.91; 2.02; 4.99; 9.84; 2.76; 2.83; 4.45; 2.49; 11.73; 4.38; 24.14; 19.89]

const obs_inc_dist = IncomeDist(inc_levels,inc_weights)
#= Initial set up of product data. Need to:
  1) import data
  2) initilize a product object for each product in data
  3) Construct markets out of products
=#

function data_setup()
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
        println("Working with market ", m)
      mkt_df = prod_char_df[(prod_char_df[:,:mkt] .== m),:] # get products in market
      mkt_vec = Liquor[]
      mkt_ps_lookup = Dict{Liquor,Tuple{Array{Float64,1},Array{Float64,1}}}() # dictionary with keys as products and tuple of prices and ffs as values.
      for row in eachrow(mkt_df)
          print(".")
        # constructing product objects
        id = row[:product]
        size_d_vec = [row[:d_s_750ML],row[:d_s_1L],row[:d_s_175L]] # get size dummy values
        s_util = 0.0 # sum(size_d_vec.*size_utils) # get group utility
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
            obs_cutoff_q = obs_cutoff_q .* 12.0
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
      println("")
    end
    return markets_array
end

export size_utils, prod_utils, nested_coefs, obs_inc_dist
export data_setup

end
