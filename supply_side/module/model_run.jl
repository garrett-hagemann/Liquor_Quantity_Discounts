include("LiquorQD.jl")
using LiquorQD, DataFrames, DataArrays

#= Set up of coefficients. Product and group utility goes in product b/c those
don't multiply out well. The remainder go in a DemandCoefs object with an associated
weight =#
group_utils = [-.9980637;.8108195;-.2365993;-.7477314] # order is gin,vod,rum,sch
prod_utils= [-3.121894;-3.121704;-4.801027;-5.014518;-4.846033;-3.435482;
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
coefs_100_plus = DemandCoefs(.028152,-.2519711, -.0019593)

weights = [.2055;.3970;.2059;.1917]

coef_array = [coefs_lt30;coefs_30_to_70;coefs_70_to_100;coefs_100_plus]

#= Initial set up of product data. Need to:
  1) import data
  2) initilize a product object for each product in data
  3) Construct markets out of products
=#

prod_char_df = readtable("E:/Box/Box Sync/Snuffles Backup/gh8728/Projects/Liquor/Data/Nielsen_Panel/analysis/product_chars.csv")
mkt_ids = Array(levels(prod_char_df[:,:mkt])) # get market ids
prod_array = Liquor[]
markets_array = Market[]
for m in mkt_ids
  mkt_df = prod_char_df[(prod_char_df[:,:mkt] .== m),:] # get products in market
  mkt_vec = Liquor[]
  for row in eachrow(mkt_df)
    id = row[:product]
    group_d_vec = [row[:d_g_gin],row[:d_g_vod],row[:d_g_rum],row[:d_g_sch]] # get producrt group dummy values
    g_util = sum(group_d_vec.*group_utils) # get group utility
    prod = Liquor(id,row[:price],row[:price],row[:imported],row[:proof],g_util,prod_utils[id])
    push!(mkt_vec,prod) # push to market
    push!(prod_array,prod) # push to total product array
  end
  date_str = mkt_df[1,:date_m]
  mkt_y = parse(Int64,date_str[1:4])
  mkt_m = parse(Int64,date_str[6:end])
  tmp_mkt = Market(mkt_vec,mkt_y,mkt_m)
  push!(markets_array,tmp_mkt)
end
