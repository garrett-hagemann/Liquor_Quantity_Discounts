capture log close
clear
set more off

log using "data_explore_licenses.txt", text replace

local files: dir "./licenses_by_county" files "*.txt"

local i = 0
foreach file in `files'{
	import delimited using "licenses_by_county/`file'", delimiters("\t") rowrange(2) varnames(2) ///
		stringcols(_all) clear
	tempfile s`i'
	save `s`i''
	local i = `i' + 1
}

local i = 0
foreach file in `files'{
	append using `s`i''
	local i = `i' + 1
}

foreach var of varlist _all {
	replace `var' = stritrim(`var')
	replace `var' = strtrim(`var')
}

/* keeping only liquor stores */

keep if type == "L"

/* removing zip+4 to just zip */
replace zip = substr(zip, 1,5) if strlen(zip) > 5

save ny_liquor_licenses_2016_06, replace

// stores by zip
preserve
	contract zip
	tab _freq
restore

// stores in NYC
preserve
	keep if county == "BRON" | county == "RICH" | county == "QUEE" | county == "KING" | county == "NEW"
	count
restore

log close

