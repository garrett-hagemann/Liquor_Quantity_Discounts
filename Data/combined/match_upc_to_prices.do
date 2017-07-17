capture log close
clear
set more off

log using "match_upc_to_prices.txt", text replace

/* This file matches the nielsen UPCs to the NYS price schedule data.
This uses the hand-checked UPC crosswalk file to do the match. The
match quality depends on the extent of the mathing. That is, as more 
products are matched by hand the quality does up. */

/* importing hand matches */
import delimited using upc_crosswalk_checked.csv
keep if checked == 1 // keeping hand checked records. Ignores multiple matches

drop similscore neg_similscore brand_string1

/* merging top upc info back in */

merge 1:m upc_id brand_string using upc_file, keep(2 3) gen(_merge_matched)

local price_sched_file = "../NY_price_schedules/old_format_months.dta"
preserve
	use `price_sched_file', clear
	drop _merge
	reshape wide disc_p disc_q disc_dollars actual_p tariff norm_p norm_tariff, i(record_num) j(option)
	tempfile ps
	save `ps'
restore

merge m:1 record_num using `ps', keep(1 3) gen(_merge_ps)

save merged_sample, replace
export delimited merged_sample.csv, replace

log close
