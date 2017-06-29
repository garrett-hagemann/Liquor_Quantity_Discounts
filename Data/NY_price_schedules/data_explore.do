clear
capture log close
set more off

log using "data_explore.txt",text replace

use old_format_months

drop _merge
drop if posting_year == 2014 // only incomplete data for 2014. No need to keep it.

replace option=option+1 // correcting options which start at 0
egen total_options = max(option), by(record_num) // calculating total number of options

// Need to reshape to one record per row for some operations. Doing that first so everything is consistent
reshape wide disc_p disc_q disc_dollars actual_p tariff norm_p norm_tariff, i(record_num) j(option)
egen product = group(size brand_name brand_label_reg_id nv_serial_number) // "product" is a brand label id and a wholesaler id
gen str_month = string(posting_month+1,"%02.0f")
gen str_year = string(posting_year)
gen str_date = str_year + "-" + str_month
gen date_n = monthly(str_date,"YM")
format date_n %tm
gen month_n = month(dofm(date_n))
gen year_n = year(dofm(date_n))

/* dropping some duplicates. Appears to be due to weak definition of a "product"
given the variables in a data set. Works mostly - but some producers appear to 
overlap numbers */
duplicates tag product date_n, gen(tag)
drop if tag > 0

reshape long // going back to long for a few things

table total_options option, contents(mean actual_p sd actual_p med actual_p) format(%9.2f)
table total_options option, contents(mean disc_q sd disc_q med disc_q) format(%9.2f)

// now same tables with easy latex integration
tabout total_options option, 

reshape wide

// distribution of number of price segments
hist total_options, discrete percent xtitle("Price Schedule Parts")
graph export option_hist.pdf, replace

tab total_options

xtset product date_n

// number of case discounts vs bottle discounts
tab discount_size_type

// price schedules per  month
tab year_n month_n
tab year_n month_n if discount_size_type == "CASE"
tab year_n month_n if discount_size_type == "BOTTLE"

// analyzing time differences

gen p1_diff = actual_p1 - L.actual_p1
gen p2_p1spread = actual_p1 - actual_p2
gen spread_diff = p2_p1spread - L.p2_p1spread

tabstat disc_q2 p2_p1spread p1_diff spread_diff if discount_size_type == "CASE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)
tabstat disc_q2 p2_p1spread p1_diff spread_diff if discount_size_type == "BOTTLE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)

gen p1_diff_pct = (actual_p1 - L.actual_p1)/L.actual_p1
gen p2_p1spread_pct = (actual_p1 - actual_p2)/actual_p1
gen spread_diff_pct = p2_p1spread_pct - L.p2_p1spread_pct

tabstat p2_p1spread_pct p1_diff_pct spread_diff_pct if discount_size_type == "CASE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)
tabstat p2_p1spread_pct p1_diff_pct spread_diff_pct if discount_size_type == "BOTTLE", stat(p1 p5 p25 p50 p75 p95 p99 mean) columns(statistics)

log close
