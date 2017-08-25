/* This file starts the data exploration on the NYS liquor price schedules.
The files are all the same format up through the middle of 2014. After which
there is a format change. The steps are:

1) get old format in and cleaned up for all months
2) get new format in and cleaned up for all months
3) combine only relevant fields (brand IDs, qty discoutns, etc.) between formats

*/

capture log close
clear all
set more off

ssc install moss // pack used for cleaning up string data. Needs to be installed

log using "data_clean_old_format.txt", text replace

/* Dealing with all years prior to 2014 */

local months "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"
local years "2009 2010 2011 2012 2013"

foreach y in `years'{
	foreach m in `months'{
		disp "Parsing `y'`m'"
		import delimited using "`y'/LR-`m'-`y'.txt", clear // will generate an error, but appears to have no effect. In stata 14, use option encoding(UTF-16). Files are encoded in UTF-16LE. But doesn't appear to contain any special characters?
		tempfile tmp`y'`m'
		save `tmp`y'`m''
	}
}

local 2014MonthsOld "Jan Feb Mar Apr"

foreach m in `2014MonthsOld'{
	disp "Parsing 2014`m'"
	import delimited using "2014/LR-`m'-2014-old-format.txt", clear // see previous UTF comment.
	tempfile tmp2014`m'
	save `tmp2014`m''
}

clear // starting with a clean dataset

foreach y in `years'{
	foreach m in `months'{
		disp "Appending `y'`m'"
		append using `tmp`y'`m''
	}
}

foreach m in `2014MonthsOld'{
	disp "Appending 2014`m'"
	append using `tmp2014`m''
}

gen record_num = _n

rename qty_discount_id discount_code

save old_format_months, replace

clear
/* importing the discount codes and then merging them in */
foreach y in `years'{
	disp "Parsing discoutn codes for `y'"
	import delimited using "`y'/Discount-codes-`y'.txt", clear // see previous UTF comment
	duplicates drop _all, force // there appear to be duplicate records
	
	tempfile discount`y'
	save `discount`y''
}

import delimited using "2014/Discount-codes-Jan-thru-Apr-2014-old-format.txt", clear // see previous UTF comment
duplicates drop _all, force // there appear to be EXACT duplicate records. Not sure why

foreach y in `years'{
	disp "Appending discount codes for `y'"
	append using `discount`y''
}

drop if posting_type == "LW"

/* again, more duplicates. Not sure why. An eye scan suggests they're the same. Perhaps there's hidden characters.
Need to be careful about this. */
duplicates drop nv_serial_number posting_year posting_month posting_type discount_code, force // droping unwated records. May want these later

tempfile AllDiscountCodes
save `AllDiscountCodes'


clear

use old_format_months

merge m:1 nv_serial_number posting_year posting_month posting_type discount_code using `AllDiscountCodes'

drop if _merge == 2 //dropping unused discount codes
/* some records appear to refer to the same product but have different prices.
Nothing in the data disctinguishes them. So, dropping for now */
duplicates drop nv_serial_number posting_year posting_month posting_type discount_code, force // droping unwated records. May want these later

save old_format_months, replace

/* Fixing sizes so they have a common format */
foreach non_char in "( ) , - . /" {
	replace size = subinstr(size, "`non_char'", "",.) // removing known non alpha-numeric characters
}
// replace size = ustrregexra(size, "[^a-zA-Z0-9.]","") // removing non-alpha numeric characters. Only works with unicode

/* cleaning up some different labels to have consistent structure. Covers 90+% of raw
data */
replace size = "1L" if size == "1.0L"
replace size = "750ML" if size == "750"
replace size = "50ML" if size == "050ML"
replace size = "1.75L" if size == "1.75"
replace size = "375ML" if size == "375"
replace size = "1L" if size == "LT"
replace size = "750ML" if size == "750BTL"
replace size = "750ML" if size == "750 BTL"
replace size = "750ML" if size == "750ML."
replace size = "1L" if size == "LITER"
replace size = "50ML" if size == "50"
replace size = "1L" if size == "1LT"
replace size = "750ML" if size == "0.75"
replace size = "200ML" if size == "200"
replace size = "750ML" if size == ".750LITERSBOTTLE"
replace size = "1L" if size == "1"
replace size = "1L" if size == "1LITER"
replace size = "1L" if size == "1LIT"
replace size = "750ML" if size == "0.75BTL"
replace size = "1LTR" if size == "1L"
replace size = "1.75L" if size == "1.75LT"
replace size = "50ML" if size == "50ML."
replace size = "375ML" if size == "375ML."
replace size = "1L" if size == "1000ML"
replace size = "1.75L" if size == "1.75LIT"
replace size = "375ML" if size == ".375LITERSBOTTLE"
replace size = "200ML" if size == ".200LITERSBOTTLE"
replace size = "1.75L" if size == "1.75ML"
replace size = "1L" if size == "1.000LITERSBOTTLE"

// catching a few more
replace size = subinstr(size, "LITER", "L",.)
replace size = subinstr(size, "LTR", "L",.)
replace size = subinstr(size, "LT", "L",.)
replace size = subinstr(size, "LIT", "L",.)

// removing spaces to standardize
replace size = subinstr(size, " ", "",.)

// one last one
replace size = "750ML" if size == "750BTL"

/* keeping only larger bottles. Most likely to be purcahsed sizes*/
keep if size == "750ML" | size == "1L" | size == "1.75L" | size == "375ML"

/* identifying the type of discount, case or bottle */
gen discount_size_type = ""
replace discount_size_type = "CASE" if regexm(discount_code, "(CS)|(CASE)|(CAS)")
replace discount_size_type = "BOTTLE" if regexm(discount_code, "(BT)|(BTL)|(BOTTLE)")

/* Identifying if discount is in percentage terms or dollar terms */
gen pct_discount = regexm(discount_code, "[\%]") // in brackets to ensure it matches %
gen dollar_discount = regexm(discount_code, "[\$]")

/* Dropping unclear discounts */
drop if pct_discount == 0 & dollar_discount == 0 // unclear discount style
drop if discount_size_type == "" // unclear if case or bottles
drop if pct_discount == 1 & dollar_discount == 1 // small number of observations with typos. Can fix later

replace discount_code = subinstr(discount_code, ")", ",",.) // normalizing some formatting
replace discount_code = subinstr(discount_code, "(", "",.) // normalizing some formatting
replace discount_code = stritrim(discount_code) // removing internal double spaces

drop if regexm(discount_code, "-") // irregular formatting for discoutn code. Maybe fix later.

/* Using the SSC package moss to pull numbers out of discount codes */

moss discount_code, match("([0-9\.]+)") regex // pulls all numbers out of discount code

/* Now converting these into useful variables. The discount descriptions always
has two paired numbres. The first is the discount amount and the second is the
required purchase quantity. So, starting with 1 it should alternate between
the discoutn amount and qty. */

foreach v of varlist _match*{
	local num : subinstr local v "_match" "" // removes "_match" form local v
	if mod(`num',2) != 0{
		local new_num = (`num' + 1) / 2
		rename `v' disc_p`new_num'
	}
	else {
		local new_num = (`num')/2
		rename `v' disc_q`new_num'
	}
}

/* Converting to numeric */
destring(disc_*), replace force // shouldn't be any non-numerics. Not sure why force is needed

drop _count _pos*

/* Restructuring into long format. This will allow graphing of the tariff,
price paid, and normalizing to see curvature. */

/* adding in a 0th option which is the normal price */
gen disc_p0 = 0
gen disc_q0 = 0

reshape long disc_p disc_q, i(record_num) j(option)

/* Dropping observations that are extra (i.e. have 2 tariff options while
the maximum observed in the data set is higher than 2) or some other problem
that resulted in a missing value for p or q */
drop if disc_p == . | disc_q == .

/* Calculating actual price paid at each point. 4 cases to deal with:
1) pct discount on cases
2) pct discount on bottles
3) dollar discount on cases
4) dollar discount on bottles

First thing is to convert all discounts into one type. So I'll convert everything
into dollar discounts. Everything will be kept in the same units as the discount
size type for now (i.e. if its a bottle discount then the tariff will be expressed
in bottle terms).

Generally there doesn't appear to be much of a discount for buying a case vs 
bottles striaght up (beyond the discount for a cases worth of bottles). In those
instances where there is a big discrepancy represent less than 1% of the data, so
we can probably assume their typos */
gen disc_dollars = .
replace disc_dollars = disc_p if dollar_discount == 1
replace disc_dollars = disc_p/100*price_per_case if discount_size_type == "CASE" & pct_discount == 1
replace disc_dollars = disc_p/100*price_per_bottle if discount_size_type == "BOTTLE" & pct_discount == 1

gen actual_p = .
replace actual_p = price_per_bottle - disc_dollars if discount_size_type == "BOTTLE"
replace actual_p = price_per_case - disc_dollars if discount_size_type == "CASE"

gen tariff = actual_p * disc_q

egen max_p = max(actual_p), by(record_num)

gen norm_p = actual_p / max_p

gen norm_tariff = norm_p*disc_q

/* some data validation. Some records have options labeled incorrectly (i.e. non-
consecutively. Generating a test and then dropping records that don't conform */

bys record_num (option): gen test = _n
gen test_tag = (option != (test-1))
egen test_total = total(test_tag), by(record_num)
drop if test_total > 0
drop test test_tag test_total

save old_format_months, replace

log close

