cd "Nielsen_Panel/analysis"

do liquor_purchases.do

do blp_sample.do

cd "../../combined"

do NYS_UPC_match.do

// now need to hand match

do match_upc_to_prices.do
