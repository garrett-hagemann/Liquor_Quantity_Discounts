clear
capture log close

log using "data_explore.txt",text replace

use old_format_months

drop _merge

// distribution of number of price segments
preserve
	sort record_num option
	bys record_num: keep if _n == _N
	replace option = option+1
	hist option, discrete percent xtitle("Price Schedule Parts")
	graph export option_hist.pdf, replace
restore

// trying to pin down how schedules change over time
reshape wide disc_p disc_q disc_dollars actual_p tariff norm_p norm_tariff, i(record_num) j(option)
egen product = group(size brand_name brand_label_reg_id nv_serial_number) // "product" is a brand label id and a wholesaler id
gen str_month = string(posting_month+1,"%02.0f")
gen str_year = string(posting_year)
gen str_date = str_year + "-" + str_month
gen date_n = monthly(str_date,"YM")
format date_n %tm

/* dropping some duplicates. Appears to be due to weak definition of a "product"
given the variables in a data set. Works mostly - but some producers appear to 
overlap numbers */
duplicates tag product date_n, gen(tag)
drop if tag > 0

xtset product date_n

gen p0_diff = actual_p0 - L.actual_p0
gen p1_p0spread = actual_p0 - actual_p1
gen spread_diff = p1_p0spread - L.p1_p0spread

tabstat disc_q1 p1_p0spread p0_diff spread_diff if discount_size_type == "CASE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)
tabstat disc_q1 p1_p0spread p0_diff spread_diff if discount_size_type == "BOTTLE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)

gen p0_diff_pct = (actual_p0 - L.actual_p0)/L.actual_p0
gen p1_p0spread_pct = (actual_p0 - actual_p1)/actual_p0
gen spread_diff_pct = p1_p0spread_pct - L.p1_p0spread_pct

tabstat p1_p0spread_pct p0_diff_pct spread_diff_pct if discount_size_type == "CASE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)
tabstat p1_p0spread_pct p0_diff_pct spread_diff_pct if discount_size_type == "BOTTLE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)

log close
