capture log close
clear
set more off
version 14

log using "ind_logit_est.txt", text replace


use "~/Liquor/Data/Nielsen_Panel/nielsen_extracts/analysis/individual_logit_sample.dta"

/* making some variables */

// Month & Year dummies
gen date_n = date(purchase_date, "YMD")
format date_n %td
gen purchase_year = year(date_n)
gen purchase_month = month(date_n)

// Liquor type dummies
gen d_gin = (product_module_descr == "GIN")
gen d_vod = (product_module_descr == "VODKA")
gen d_rum = (product_module_descr == "RUM")
gen d_sch = (product_module_descr == "SCOTCH")
gen d_brb = regexm(product_module_descr, "BOURBON")
gen d_whs = regexm(product_module_descr, "WHISKEY")

// generating choice variable. Should be 0 for non-liquor trips
egen choice = group(upc upc_ver_uc) if liquor_trip == 1 // may want to just group by upc instead of versions
replace choice = 0 if liquor_trip == 0  // setting choice to 0 for non-liquor trips

gen price = total_price_paid / quantity
replace price = 0 if liquor_trip == 0 // setting price to 0 for non-liquor trips

keep purchase_year purchase_month upc upc_ver_uc d_* choice price

// approximating the top number of choices
egen _freq = count(choice), by(choice)

egen freq_rank = rank(_freq), field

egen rank_rank = group(freq_rank)


//defining constraints
local rank_thresh = 55
local rhs "price d_gin d_vod d_rum d_sch d_brb d_whs"
levelsof choice if rank_rank < `rank_thresh' & rank_rank > 1, local(choices) // levels of choice except the 0

local counter = 1
local first : word 1 of `choices'
local rest : list choices - first
foreach c of local rest{
	constraint `counter++' [`first' = `c']: `rhs'
} 

preserve
	contract choice `rhs' rank_rank, freq(_est_freq)
	count
	mlogit choice `rhs' [fweight = _est_freq] if rank_rank < `rank_thresh', constraints(1/`counter')
restore 

/* unfortunately this fails to converge. Stata insists there is a discontinuous region. This is likely because
while this should match the underlying utility function, Stata probably doesn't like to be forced to estimate
this way. If you consider the results in logit_vs_reg_test.do it looks like there may be a linear way to do it.
This would be equivalent to the Berry-Logit proceedure - just IV regression on the adjusted market shares. 
The advantage of the logit-style approach is the ability to include individual level regressors (which would 
then have choice specific effects). However, including instruments is more difficult in the logit style approach. 
Perhaps there is a GMM style approach to the individual level logit that could include instruments. Regardless,
it appears to be key to reduce the size of the dataset and then frequency weight. The dataset goes from
3.2 million to like 16k. This is especially helpful given that most of the observations are non-liquor trips
which are the outside option. */

log close
