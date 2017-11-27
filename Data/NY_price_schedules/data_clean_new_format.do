/* This file cleans the new format of the NY liquor price schedules.
This constitutes the second half of 2014 and all of 2015. Once cleaned,
subsets of the data can be stacked with the old format with a little renaming.*/

capture log close
version 13
set more off
clear all
log using data_clean_new_format.txt, text replace

local months "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"


foreach m in `months'{
	disp "Parsing 2015`m'"
	import delimited using "2015/LR-`m'-2015.txt", clear stringcols(_all)
	tempfile tmp2015`m'
	save `tmp2015`m''
}

local 2014MonthsNew "May Jun Jul Aug Sep Oct Nov Dec"

foreach m in `2014MonthsNew'{
	disp "Parsing 2014`m'"
	import delimited using "2014/LR-`m'-2014-new-format.txt", clear stringcols(_all)
	tempfile tmp2014`m'
	save `tmp2014`m''
}

clear // starting with a clean dataset

foreach m in `months'{
	disp "Appending 2015`m'"
	append using `tmp2015`m''
}


foreach m in `2014MonthsNew'{
	disp "Appending 2014`m'"
	append using `tmp2014`m''
}

gen record_num = _n
destring *_price amount_* qty_*, replace // converting to numeric
drop if amount_2 == . // set of 16 records appears to have garbled input. Can fix later if need be.

destring *_price amount_* qty_*, replace // trying again

destring post_month post_year, replace

// correcting item sizes
replace item_size = upper(item_size) // remove lower case
replace item_size = "750ML" if item_size == "750"
replace item_size = "750ML" if item_size == "750.00"
replace item_size = "1L" if item_size == "1.00"
replace item_size = "1.75L" if item_size == "1.75"
replace item_size = "1L" if item_size == "1"
replace item_size = "375ML" if item_size == "375"
replace item_size = "375ML" if item_size == "375.00"
replace item_size = "50ML" if item_size == "50"
replace item_size = "50ML" if item_size == "50.00"
replace item_size = "200ML" if item_size == "200"
replace item_size = "200ML" if item_size == "200.00"
replace item_size = "750ML" if item_size == "0.75"
replace item_size = "750ML" if item_size == "0.75L"
replace item_size = "1L" if item_size == "1LT"
replace item_size = "100ML" if item_size == "100"
replace item_size = "1L" if item_size == "1.0L"
replace item_size = "1.75L" if item_size == "1.75LT"
replace item_size = "100ML" if item_size == "100.00"
replace item_size = "1L" if item_size == "1.0"
replace item_size = "1L" if item_size == "1LIT"
replace item_size = "1.75L" if item_size == "1.75LIT"
replace item_size = "1L" if item_size == "1000ML"
replace item_size = "1L" if item_size == "1000"

keep if item_size == "750ML" | item_size == "1L" | item_size == "1.75L" | item_size == "375ML" // keeping only standard sizes

// replace quantity == . if no disocunt is offered

foreach v of varlist qty_1 qty_2 qty_3 qty_4 qty_5 qty_6 qty_7 qty_8 qty_9 qty_10 {
	replace `v' = . if `v' == 0
}

foreach v of varlist unit_1 unit_2 unit_3 unit_4 unit_5 unit_6 unit_7 unit_8 unit_9 unit_10 {
	replace `v' = "C" if !(`v' == "B" | `v' == "")
}

/* Generating "price paid" at each quantity since discounts are given by the 
amount TAKEN OFF of the price per case */
rename qty_discid disc_descr
reshape long qty_ unit_ amount_ type_, i(record_num) j(option) string

gen actual_p = .
replace actual_p = case_price*(1-amount) if unit == "C" & type == "%" & qty != .
replace actual_p = bot_price*(1-amount) if unit == "B" & type == "%" & qty != .
replace actual_p = case_price-amount if unit == "C" & type == "$" & qty != .
replace actual_p = bot_price-amount if unit == "B" & type == "$" & qty != .

reshape wide qty_ unit_ amount_ actual_p type, i(record_num) j(option) string

keep if post_year == 2015

save new_format_months, replace
