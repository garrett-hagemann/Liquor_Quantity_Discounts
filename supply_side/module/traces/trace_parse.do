capture log close
clear
set more off

log using trace_parse.txt, text replace

/* generating list of files */
local traces : dir "." files "*.csv"

local i = 0
foreach t of local traces {
	import delim "`t'", clear
	rename v1 c
	rename v2 b 
	rename v3 Q
	
	foreach var of varlist c b Q {
		egen `var'_mode = mode(`var')
		egen `var'_med = median(`var')
		egen `var'_mean = mean(`var')
		egen `var'_sd = sd(`var')
	}
	
	keep *_* // keeping summary stats
	keep in 1 // keeping just one obs

	gen year = regexs(1) if regexm("`t'", "trace_([0-9]+)_([0-9]+)_([0-9]+).csv")
	gen month = regexs(2) if regexm("`t'", "trace_([0-9]+)_([0-9]+)_([0-9]+).csv")
	gen prod = regexs(3) if regexm("`t'", "trace_([0-9]+)_([0-9]+)_([0-9]+).csv")
	
	destring _all, replace
	
	local save_str : subinstr local t ".csv" "", all
	
	tempfile tmp`save_str'
	save `tmp`save_str''
}

clear

foreach t of local traces {
	local save_str : subinstr local t ".csv" "", all
	append using `tmp`save_str''
}

kdensity c_mode
graph export c_density.pdf, replace

kdensity b_mode
graph export b_density.pdf, replace


destring _all, replace

save trace_sum, replace

log close
