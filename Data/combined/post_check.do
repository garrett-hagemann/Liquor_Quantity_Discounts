capture log close
clear
set more off

log using "post_check.txt", text replace

/* This file matches the nielsen UPCs to the NYS price schedule data.
This uses the hand-checked UPC crosswalk file to do the match. The
match quality depends on the extent of the mathing. That is, as more 
products are matched by hand the quality does up. */

/* importing hand matches */
import delimited using upc_crosswalk_checked.csv
keep if checked == 1 // keeping hand checked records. Ignores multiple matches

drop similscore neg_similscore brand_string1

// fixing bad matching. Just need to match on brand_string and month year
format upc_id %12.0f
gen purchase_month = regexs(1) if regexm(string(upc_id,"%12.0f"),"^([0-9]+)(2[0-1)[0-9][0-9])([0-9]+)$")
destring purchase_month, replace
gen purchase_year = substr(string(upc_id,"%12.0f"),2,4) if purchase_month < 10
replace purchase_year = substr(string(upc_id,"%12.0f"),3,4) if purchase_month >= 10
destring purchase_year, replace

/* merging top upc info back in */

merge 1:m purchase_month purchase_year brand_string using upc_file, keep(2 3) gen(_merge_matched)
drop if product == 0

local price_sched_file = "../NY_price_schedules/old_format_months.dta"
preserve
	use `price_sched_file', clear
	drop _merge
	reshape wide disc_p disc_q disc_dollars actual_p tariff norm_p norm_tariff, i(record_num) j(option)
	tempfile ps
	save `ps'
restore

merge m:1 record_num using `ps', keep(1 3) gen(_merge_ps)

gen type = ""
replace type = "oo" if product == 0
replace type = "gin" if d_g_gin
replace type = "vod" if d_g_vod
replace type = "rum" if d_g_rum
replace type = "sch" if d_g_sch
replace type = "brb" if d_g_brb
replace type = "whs" if d_g_whs
replace type = "teq" if d_g_teq
replace type = "otr" if d_g_otr


/*
/* hack because match was done with only 101 products*/
local J = 250
replace product = (`J'+2) if product == (`J'+1) & d_g_vod == 1
replace product = (`J'+3) if product == (`J'+1) & d_g_sch == 1
replace product = (`J'+4) if product == (`J'+1) & d_g_brb == 1
replace product = (`J'+5) if product == (`J'+1) & d_g_whs == 1
replace product = (`J'+6) if product == (`J'+1) & d_g_teq == 1
replace product = (`J'+7) if product == (`J'+1) & d_g_otr == 1
replace product = (`J'+8) if product == (`J'+1) & d_g_rum == 1
*/

/* taking care of duplicate product records in a month */
//duplicates drop product date_m if product <= 100, force

sort date_m product

egen mkt = group(date_m) // doesnt matter if it doesn't match mkt definition elsewhere. No mkt specific variables

//bys date_m product (_merge_ps): keep if _n == _N // only need one record of product per mkt. Not per case. Keeps matched records

save merged_sample, replace
export delimited merged_sample.csv, replace

log close
