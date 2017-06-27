capture log close
clear
set more off


log using "blp_sample.txt", text replace

/* This file makes the relevant variables for using BLP estimation
proceedure found in demand_est/BLP-basic. Due to the relatively low
market shares we're using the state wide purcahses. This makes each
"market" a month-year combo, coinciding with the price schedule setting
period.

Currently focusing on just a few UPCs to get the code working. Can
continue to add more products. Likewise, want to include price sched
information as instruments for price.


TODO:
	Fix HACK
*/

use ny_liquor_purchases_all_years

/* Defining markets */
egen mkt = group(purchase_year purchase_quarter fips_county_code)

/*
/* Defining Market Size  according to number of purchases */
preserve
	clear
	foreach year in 2009 2011 2012 2013 2014 {
		append using ny_trips_`year'
	}
	gen date_n = date(purchase_date, "YMD")
	gen purchase_year = year(date_n)
	gen purchase_month = month(date_n)
	
	contract purchase_year purchase_month panelist_zip_code, freq(M)

	tempfile mkt_size
	save `mkt_size'
restore

merge m:1 purchase_year purchase_month panelist_zip_code using `mkt_size', gen(_merge_mkt) keep(3)
*/

/* generating county dummies */
gen d_c_bronx = (fips_county_descr=="BRONX")
gen d_c_kings = (fips_county_descr=="KINGS")
gen d_c_ny = (fips_county_descr=="NEW YORK")
gen d_c_queens = (fips_county_descr=="QUEENS")
gen d_c_richmond = (fips_county_descr=="RICHMOND")

/* Defining market size based on Census. 2010 population in county >21 years. The
implication of this, and the quarterly frequency, is that people buy a bottle of
liquor once a quarter on average */

gen M = .
replace M = 944003 if fips_county_descr=="BRONX"
replace M = 1804542 if fips_county_descr=="KINGS"
replace M = 1289661 if fips_county_descr=="NEW YORK"
replace M = 1681493 if fips_county_descr=="QUEENS"
replace M = 340148 if fips_county_descr=="RICHMOND"

replace M = M*2 // changing to buying 2 bottles a quarter on average.
// TODO: Should figure out scaling from Nielsen purchase frequency

/* Generating additional variables */

egen product = group(upc) // maybe use UPC ver as well? may not matter much
gen price = total_price_paid / quantity

/* Generating quantity weighted average price. Given that
total_price_paid = unit_price*qty we can just sum that for numerator and take
total quantity as the denominator */

egen total_tpp = total(total_price_paid), by(product mkt)
egen total_q = total(quantity), by(product mkt)
gen avg_price = total_tpp / total_q // quantity weighted average price

egen total_upc_purchases = total(quantity), by(product)
drop if total_upc_purchases < 10 // dropping products not observed frequently

egen share = total((quantity*projection_factor)), by(product mkt) // total purchases of a UPC by mkt
replace share = share / M // converting total purchases to shares

/* Need to get rid of duplicate purchases of the same product in a market at
this point. They were needed for the share calculation but should be removed
at this point so as not to double count shares. Also removing erroneous records
here. 0 prices, unwanted categories, etc. */

bys mkt product: keep if _n == 1
drop if product_module_descr == "COOLERS - REMAINING" // dropping from analysis
drop if product_module_descr == "ALCOHOLIC COCKTAILS" // dropping from analysis
drop if avg_price == 0 // not sure what these are. Likely data entry errors
keep if size1_amount == 1 | size1_amount == 1.75 | size1_amount == 375 ///
	| size1_amount == 750 // keeping 1L, 1.75L, 375ML, 750ML sizes

egen inside_share = total(share), by(mkt)
gen outside_share = 1 - inside_share

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



order mkt M purchase_year purchase_month purchase_quarter product price share inside_share ///
	outside_share d_* proof
	

save blp_sample, replace

export delimited using "blp_sample.csv", delim(",") quote replace

log close
