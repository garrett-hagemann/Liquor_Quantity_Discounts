capture log close
clear
set more off

log using cf_parse.txt, text replace

/* generating list of files */
local cfs : dir "cf/" files "*.csv"

foreach cf of local cfs {
	disp "`cf'"
	import delim "cf/`cf'", clear
	
	local save_str : subinstr local cf ".csv" "", all
	
	tempfile tmp`save_str'
	save `tmp`save_str''
}

disp "All files parsed. Combining."

clear

foreach cf of local cfs {
	local save_str : subinstr local cf ".csv" "", all
	append using `tmp`save_str''
}


merge 1:1 month year prod using trace_sum
drop _merge
rename prod product
rename month purchase_month
rename year purchase_year

merge 1:1 product purchase_month purchase_year using "../../Data/combined/merged_sample", gen(_merge_cf)

save cf_sum, replace

log close
