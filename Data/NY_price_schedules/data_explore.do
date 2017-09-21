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
egen product = group(size brand_name product_item_name nv_serial_number) // "product" is a brand label id and a wholesaler id
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

// for cases
table total_options option if discount_size_type == "CASE", contents(mean actual_p sd actual_p med actual_p freq) format(%9.2f)
table total_options option if discount_size_type == "CASE", contents(mean disc_q sd disc_q med disc_q freq) format(%9.2f)

// now same tables with easy latex integration
//tabout total_options option using tab_options.txt, replace


reshape wide

// distribution of number of price segments
hist total_options if discount_size_type == "CASE", discrete percent xtitle("Price Schedule Parts")
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

preserve
	egen product2 = group(size brand_name product_item_name)
	contract product2  posting_month posting_year
	tab _freq
restore


/* based on my own analysis. Has codes that don't appera in file formats */
gen spirit_type = ""
replace spirit_type = "vod" if beverage_type == "A"
replace spirit_type = "sch/whs/brb" if beverage_type == "B"
replace spirit_type = "gin" if beverage_type == "C"
replace spirit_type = "brandy" if beverage_type == "D"
replace spirit_type = "rum" if beverage_type == "E"
replace spirit_type = "liqueur" if beverage_type == "F"
replace spirit_type = "cocktail" if beverage_type == "G"
replace spirit_type = "otr" if beverage_type == "K"
replace spirit_type = "teq" if beverage_type == "M"
replace spirit_type = "mislabeled" if beverage_type == "O"

encode spirit_type, gen(spirit_type_n)
encode size, gen(size_n)

drop if spirit_type == "mislabeled" | spirit_type == "cocktail" | spirit_type == "otr"

poisson total_options actual_p1 disc_q2  proof_alcohol_cont i.spirit_type_n i.size_n i.month_n if discount_size_type == "CASE", cformat(%3.2f) pformat(%3.2f)

reshape long
reg actual_p i.option i.total_options proof_alcohol_cont i.spirit_type_n i.size_n i.month_n if discount_size_type == "CASE" & total_options <= 8 & total_options > 1


log close
