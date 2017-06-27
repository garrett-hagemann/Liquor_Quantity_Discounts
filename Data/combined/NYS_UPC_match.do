/* The goal of this do file is to match the brands observed in the NYS price
schedule data to the UPC data in the Nielsen data. Unfortunately, the strings
describing the brands are wildly different. The first step is to use some fuzzy
string matching algorithms availble in SSC. After that, any clean up will
be done by hand in this file. Any hand matching done will have a corresponding
key file to shows the linkages.

The fuzzy match uses --matchit-- which is available from SSC. It depends on
--freqindex-- which is also on SSC. Ensure that those are installed first */

capture log close
clear
version 14
set more off

log using "NYS_UPC_match.txt", text replace

/* To improve the quality of the match we're going to do this by year,month
and product size. There may be some duplicate listings in the price schedules
by these categories, due to different wholesalers. However, a quick look
shows that the difference in the prices is pretty minimal. This is likely
something that should be fixed */

local years = "2009 2010 2012 2013 2014"
local months = "1 2 3 4 5 6 7 8 9 10 11 12"
local sizes = "750ML 1L 1.75L"

/* first use the old format NYS data */
use "../NY_price_schedules/old_format_months"

drop _merge // need to drop for reshape
// reshaping wide so record num is unique key
reshape wide disc_p disc_q disc_dollars actual_p tariff norm_p norm_tariff, i(record_num) j(option)

// making one long string that is to be matched. Contains brand, product, and size
gen brand_string = brand_name + ";" + product_item_name + ";" + size

foreach year in `years'{
	foreach month in `months'{
		foreach size in `sizes'{
			if "`size'" == "1.75L"{
				tempfile ps_y`year'_m`month'_s175L
				preserve
					keep if (posting_month == `month' - 1) ///
						& (posting_year == `year') ///
						& (size == "`size'")
					
					save `ps_y`year'_m`month'_s175L' 
				restore
				disp "y`year',m`month',s`size'"
			}
			else {
				tempfile ps_y`year'_m`month'_s`size'
				preserve
					keep if (posting_month == `month' - 1) ///
						& (posting_year == `year') ///
						& (size == "`size'")
					
					save `ps_y`year'_m`month'_s`size''
				restore
				disp "y`year',m`month',s`size'"
			}
		}
	}
}

clear

use "../Nielsen_Panel/analysis/ny_liquor_purchases_all_years"

/* Now we need to separate purchases by month, year, and size to allign
with the price schedule data. Then we can apply the matching algorithm */

foreach year in `years'{
	foreach month in `months'{
		foreach size in `sizes'{
			if "`size'" == "1.75L"{
				tempfile res_y`year'_m`month'_s175L
				tempfile raw_y`year'_m`month'_s175L
				preserve
					keep if (purchase_month == `month') ///
					& (purchase_year == `year') ///
					& (size1_amount == 1.75)	

					matchit upc_id brand_string using `ps_y`year'_m`month'_s175L', idusing(record_num) txtusing(brand_string) override sim(ngram, 3)

					gen neg_similscore = -similscore

					sort upc_id neg_similscore // sorting by each UPC and then by match rank within

					save  `raw_y`year'_m`month'_s175L' // useful to fix hand matches

					by upc_id: keep if _n == 3 // keeping top 3 matches 

					save `res_y`year'_m`month'_s175L'
					
				restore
			}

			if "`size'" == "1L"{
				tempfile res_y`year'_m`month'_s1L
				tempfile raw_y`year'_m`month'_s1L
				preserve
					keep if (purchase_month == `month') ///
					& (purchase_year == `year') ///
					& (size1_amount == 1)	

					matchit upc_id brand_string using `ps_y`year'_m`month'_s1L', idusing(record_num) txtusing(brand_string) override sim(ngram, 3)

					gen neg_similscore = -similscore

					sort upc_id neg_similscore // sorting by each UPC and then by match rank within

					save  `raw_y`year'_m`month'_s1L' // useful to fix hand matches

					by upc_id: keep if _n == 3 // keeping top 3 matches 

					save `res_y`year'_m`month'_s1L'
				restore
			}

			if "`size'" == "750ML"{
				tempfile res_y`year'_m`month'_s750ML
				tempfile raw_y`year'_m`month'_s750ML
				preserve
					keep if (purchase_month == `month') ///
					& (purchase_year == `year') ///
					& (size1_amount == 750)	

					matchit upc_id brand_string using `ps_y`year'_m`month'_s750ML', idusing(record_num) txtusing(brand_string) override sim(ngram, 3)

					gen neg_similscore = -similscore

					sort upc_id neg_similscore // sorting by each UPC and then by match rank within

					save  `raw_y`year'_m`month'_s750ML' // useful to fix hand matches

					by upc_id: keep if _n <= 3 // keeping top 3 matches 

					save `res_y`year'_m`month'_s750ML'
				restore
			}
		}
	}
}

clear

local sizes = "175L 1L 750ML"

foreach year in `years'{
	foreach month in `months'{
		foreach size in `sizes'{
			append using `res_y`year'_m`month'_s`size''
		}
	}
}			

save upc_crosswalk, replace
export delim using "upc_crosswalk.csv", replace delim(",") quote

clear

foreach year in `years'{
	foreach month in `months'{
		foreach size in `sizes'{
			append using `raw_y`year'_m`month'_s`size''
		}
	}
}			

save upc_match_raw, replace

log close
