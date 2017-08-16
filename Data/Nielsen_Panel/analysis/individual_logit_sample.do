capture log close
clear
set more off

log using "individual_logit_sample.txt", text replace

local rebuild_data = 1 // set to 1 to rebuild data from scratch. Very slow
local J = 250 // number of products to get

set seed 317596364  // seed taken from random.org


local years = "2009 2010 2011 2012 2013 2014" // "2008 2009 2010 2011 2012 2013 2014"
if `rebuild_data' { 
/* Importing UPCs to keep just the liquor ones */
	import delimited using "../nielsen_extracts/HMS/Master_Files/Latest/products.tsv", stringcols(1)

	/* Keeping only Liquor records. CAn narrow it down to Alcoholic beverage first,
	but no need. Only other Alocohol is coded as Beer or Wine. There may be some
	marginal cases that differ between the NYS data and the Nielsen Data. Right now, 
	just working with the Liquor records for simplicity */
	keep if product_group_desc == "LIQUOR"  // corresponding product_group_code is 5002

	/* Uncomment the following line to remove multi-packs. This might make matching
	brand strings easier between the NYS and Nielsen datasets. The price schedules do
	contain multipacks, but it seems more complicated to model. The multipacks represent
	a form of retail quantity discounts. The focus here is on wholesale quantity discoutns

	drop if multi > 1 */

	save liquor_upcs, replace

	/* generating varibles that we want */
	//dropping unwanted types of products
	drop if product_module_descr == "COOLERS - REMAINING"
	drop if product_module_descr == "ALOCOHOLIC COCKTAILS"

	// Liquor type dummies
	gen d_g_gin = (product_module_descr == "GIN")
	gen d_g_vod = (product_module_descr == "VODKA")
	gen d_g_rum = (product_module_descr == "RUM")
	gen d_g_sch = (product_module_descr == "SCOTCH")
	gen d_g_brb = regexm(product_module_descr, "BOURBON")
	gen d_g_whs = regexm(product_module_descr, "WHISKEY")
	gen d_g_teq = (product_module_descr == "TEQUILA") 
	gen d_g_otr = (d_g_gin == 0 & d_g_vod == 0 & d_g_rum == 0 & d_g_sch == 0 & d_g_brb == 0 & d_g_whs == 0 & d_g_teq == 0)

	// Size dummies
	gen d_s_375ML = (size1_amount == 375)
	gen d_s_750ML = (size1_amount == 750)
	gen d_s_1L = (size1_amount == 1)
	gen d_s_175L = (size1_amount == 1.75) 

	/* Generating proof */
	gen proof = regexs(1) if regexm(upc_descr , "([0-9]?[0-9][0-9][.]?[0-9]?[0-9]?)P")
	destring proof, replace

	/* imported flag */
	gen imported = regexm(upc_descr , " IM ")
	replace imported = 1 if regexm(upc_descr, " PRD ") & d_g_rum == 1 // PRD is Puerto Rico Distilled
	replace imported = 1 if d_g_sch == 1 // All scotch is imported
	replace imported = 1 if d_g_teq == 1 // all tequila is imported

	// imputing proof for products which the regex didn't work
	reg proof d_g_* imported, nocons
	predict proof_hat, xb
	replace proof = proof_hat if proof == .

	tempfile small_upc_file
	save `small_upc_file'

	foreach year in `years'{
		
		use ny_trips_`year', clear
		gen NYC = (fips_county_descr == "BRONX" | fips_county_descr == "NEW YORK" | fips_county_descr == "KINGS" | fips_county_descr == "QUEENS" | fips_county_descr == "RICHMOND")
		//keep if NYC == 1
		keep trip_code_uc NYC
		tempfile tmp_trips
		save `tmp_trips', replace
				
		import delimited using "../nielsen_extracts/HMS/`year'/Annual_Files/purchases_`year'.tsv", clear colrange(1:5) stringcol(2)
		
		/* first merging in just trip codes to reduce data set size */
		
		merge m:1 trip_code_uc using `tmp_trips', gen(_merge_trips) keep(3)
		drop _merge_trips
		merge m:1 trip_code_uc using ny_trips_`year', gen(_merge_trips) keep(3)
		
		// Identifying liquor purchases
		merge m:1 upc upc_ver_uc using `small_upc_file', gen(_merge_upcs) keep(1 3) // keeps all purchases, drops unused liquor UPCs
		
		//Wine is outside option. Identifying wine purchases
		merge m:1 upc upc_ver_uc using wine_upcs, gen(_merge_wine) keep(1 3)
		
		// Keeping just one purchase for non-liquor trips
		egen max_merge = max(_merge_upcs), by(trip_code_uc) // largest merge code. Should be 3
		gen liquor_trip = (max_merge == 3) // trips where one of the UPCs is a liquor UPC
		drop if liquor_trip == 1 & _merge_upcs == 1 // dropping non-liquor records for liquor trips
		
		egen max_wine = max(_merge_wine), by(trip_code_uc)
		gen wine_trip = (max_wine == 3)
		drop if wine_trip == 0 & liquor_trip == 0 // dropping non wine non liquor trips
		
		bys trip_code_uc: gen trip_order = _n // ordering purchases within trips. Doesn't matter if sort is stable. Just need one record for each non-liquor trip
		drop if liquor_trip == 0 & trip_order > 1 // dropping additional purchases from non-liquor trips
		
		tempfile obs`year'
		save `obs`year''
	}

	clear

	foreach year in `years'{
		append using `obs`year''
	}

	
	gen date_d = date(purchase_date, "YMD") // daily date
	format date_d %td

	gen date_m = mofd(date_d) // monthly date
	format date_m %tm
	
	gen year = year(date_d)
	//keep if year == 2011

	gen price = total_price_paid / quantity

	gen white = (race == 1)
	gen oldest_age = max(male_head_age, female_head_age)

	gen choice = 1 // all these records are actual choices
	keep if d_s_375ML | d_s_750ML | d_s_1L | d_s_175L // keeping only popular size bottles

	// getting top J UPCs
	preserve
		drop if liquor_trip == 0 // keeping only liquor purchases
		contract upc, freq(upc_purchases)
		gen rev = -upc_purchases
		sort rev
		keep in 1/`J' // keeping top J products
		gen product = _n
		save top_upc_individual, replace
	restore

	merge m:1 upc using top_upc_individual, gen(_merge_top) assert(1 3)

	replace product = 0 if _merge_top == 1 & liquor_trip == 0 // outside option
	replace product = (`J'+1) if _merge_top == 1 & liquor_trip == 1 // all other products
	replace product = (`J'+2) if product == (`J'+1) & d_g_vod == 1
	replace product = (`J'+3) if product == (`J'+1) & d_g_sch == 1
	replace product = (`J'+4) if product == (`J'+1) & d_g_brb == 1
	replace product = (`J'+5) if product == (`J'+1) & d_g_whs == 1
	replace product = (`J'+6) if product == (`J'+1) & d_g_teq == 1
	replace product = (`J'+7) if product == (`J'+1) & d_g_otr == 1
	replace product = (`J'+8) if product == (`J'+1) & d_g_rum == 1

	gen case = _n // Each "purchase" is a case. Now need to construct alternatives
	egen brand_string = concat(brand_descr upc_descr size1_amount size1_units), punct(";") // used to match to price schedules
	egen mkt = group(date_m)
	save prod_chars_individual, replace
}

clear
use prod_chars_individual
/* outputing csv of product characteristics */
preserve
	drop if product == 0
	export delim using "product_chars.csv", replace
restore

/* Generating population & sample weights then sampling fewer of the outside
option purchases. This should substantially reduce the number of observations
in the data set. Note that weights apply to products not people. */

count
disp "Original case count: " r(N)
egen pop_weight = count(case), by(product)
replace pop_weight = pop_weight / _N // converting to percents

// sampling non-liquor purchases
/*
count if product > `J'
local liq_N = r(N)
sample `liq_N' if product == 0, count // sampling as many outside option purchases as non-top 100 product purchases
*/
egen sample_weight = count(case), by(product)
replace sample_weight = sample_weight / _N // converting to percents

// end of sampling code

keep case product choice price d_*_* proof imported date_m household_income white oldest_age NYC pop_weight sample_weight brand_string
egen avg_price = mean(price), by(product date_m)

/* replacing imported, proof, and size dummies with avg value for aggregate products */
foreach var of varlist d_s_* proof imported {
	egen mean_`var' = mean(`var'), by(product)
	replace `var' = mean_`var' if product > `J'
	drop mean_`var'
}

tempfile all_choices
save `all_choices'

fillin product case // generate 1 product record per choice situation
replace choice = 0 if choice == . // all new records are non-choices

// setting product 0 attributes to 0
foreach var of varlist d_*_* proof import price avg_price {
	replace `var' = 0 if product == 0 
}

// filling in case variables
foreach var of varlist household_income white NYC date_m {
	bys case (`var'): replace `var' = `var'[1] if missing(`var')
}


// filling in product variables except avg_price
/* Uses means for conglomerate products and actual value for single products */
foreach var of varlist d_g_* d_s_* proof imported {
	bys product (`var'): replace `var' = `var'[1]
}

// filling in brand string to use in price sched matching
bys product (brand_string): replace brand_string = brand_string[_N] // empty strings sorted first

// fixing average price when not missing
bys product date_m (avg_price): replace avg_price = avg_price[1]

// interpolating missing avg_price
preserve
	bys product date_m: keep if _n == 1
	sort product date_m
	xtset product date_m
	gen avg_price_inter = avg_price
	// filling forward
	bys product (date_m): replace avg_price_inter = avg_price_inter[_n-1] if avg_price_inter == .
	// filling backward to handle missings at begining of data
	gen rev = -date_m
	bys product (rev): replace avg_price_inter = avg_price_inter[_n-1] if avg_price_inter == .
	sort product date_m
	keep product date_m avg_price_inter
	tempfile ipol
	save `ipol'
restore

merge m:1 product date_m using `ipol'

/*
levelsof date_m, local(months)
foreach m of local months {
	keep if date_m == `m'
	tempfile m_choices
	tempfile m_alts
	
	save `m_choices' // original choices
	
	fillin case product
	
	replace choice = 0 if _fillin == 1
	
	keep if choice == 0 // alternatives
	save `m_alts'
	
	use `m_choices', clear
	
	sample 1, count by(product) // randomly select one product record each month
	
	drop case household_income white NYC // droping case specific vars
	
	merge 1:m product using `m_alts', update replace gen(_merge_alts)
	
	append using `m_choices'
	
	sort case product
	
	tempfile `m'data
	save ``m'data'
	clear
	use `all_choices'
}

clear
foreach m of local months{
	append using ``m'data'
}
*/



/* sampling alternatives within cases to further reduce data set size. Should be OK for Logit.
Always keeping outside option and choice */
//sample 20 if (choice == 0) & (product > 0), count by(case)


// Replacing household income where missing within a case.


compress // making dataset smaller

//drop if date_m < tm(2008m1)

save individual_logit_sample, replace


log close	

