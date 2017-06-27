capture log close
clear
set more off
set seed 23406
set emptycells drop
set matsize 800

log using ind_logit_test.txt, text replace

// using 2008 for now
use ny_liquor_purchases_2008
sample 1000, count // n individuals. up to n alternatives (alternatives are sampled again below)
egen product = group(upc) // maybe use UPC ver as well? may not matter much
gen price = total_price_paid / quantity
gen purchase_id = _n
egen case = group(household_code trip_code_uc purchase_id)

// generating product characteristic file
preserve
	keep upc-_merge_upcs product price
	tempfile prod_chars
	save `prod_chars'
restore

// generating household char file
preserve
	keep household_code-_merge_trips case
	tempfile hh_chars
	save `hh_chars'
restore

preserve 
	clear
	use `hh_chars'
	cross using `prod_chars'
	sample 10, count by(case)
	tempfile alts
	save `alts'
restore

gen choice = 1

append using `alts'

replace choice = 0 if choice == .

/* dealing with repeated alternatives. This has to do with the sampling. ideally
we'd ignore the choice and sample from among the rest, but that's hard */
duplicates tag case product, gen(tag)

drop if tag > 0 // ultimately should be able to clean up these records, but lose them for now

gen cons = 1

// Liquor type dummies
gen d_gin = (product_module_descr == "GIN")
gen d_vod = (product_module_descr == "VODKA")
gen d_rum = (product_module_descr == "RUM")
gen d_sch = (product_module_descr == "SCOTCH")
gen d_brb = regexm(product_module_descr, "BOURBON")
gen d_whs = regexm(product_module_descr, "WHISKEY")
gen d_teq = (product_module_descr == "TEQUILA") 
gen d_otr = (d_gin == 0 & d_vod == 0 & d_rum == 0 & d_sch == 0 & d_brb == 0 & d_whs == 0 & d_teq == 0)

// Size dummies
gen d_s_350ML = (size1_amount == 350)
gen d_s_750ML = (size1_amount == 750)
gen d_s_1L = (size1_amount == 1)
gen d_s_175L = (size1_amount == 1.75) 

/* Generating proof */
gen proof = regexs(1) if regexm(upc_descr , "([0-9]?[0-9][0-9][.]?[0-9]?[0-9]?)P")
destring proof, replace

/* imported flag */
gen imported = regexm(upc_descr , " IM ")
replace imported = 1 if regexm(upc_descr, " PRD ") & d_rum == 1 // PRD is Puerto Rico Distilled
replace imported = 1 if d_sch == 1 // All scotch is imported
replace imported = 1 if d_teq == 1 // all tequila is imported

gen py = price/household_income
gen iy = imported*household_income

asclogit choice py, case(case) alternatives(product) ///
	collinear



