capture log close
clear
set more off
set matsize 800
set seed 1987

log using asclogit_test.txt, text replace

use individual_logit_sample.dta

// generating dummies for income levels
gen byte inc_lt30 = (household_income <= 15)
gen byte inc_30_to_70 = (household_income > 15 & household_income <= 23)
gen byte inc_70_to_100 = (household_income == 26)
gen byte inc_100_plus = (household_income > 26)

gen p_x_inc_lt30 = price*inc_lt30
gen p_x_inc_30_to_70 = price*inc_30_to_70
gen p_x_inc_70_to_100 = price*inc_70_to_100
gen p_x_inc_100_plus = price*inc_100_plus

gen neg_p_x_inc_lt30 = -p_x_inc_lt30 - 1e-2
gen neg_p_x_inc_30_to_70 = -p_x_inc_30_to_70 - 1e-2
gen neg_p_x_inc_70_to_100= -p_x_inc_70_to_100 - 1e-2
gen neg_p_x_inc_100_plus = -p_x_inc_100_plus - 1e-2

gen byte month = month(dofm(date_m))
gen int year = year(dofm(date_m))

gen byte d_m_1 = (month == 1)
gen byte d_m_2 = (month == 2)
gen byte d_m_3 = (month == 3)
gen byte d_m_4 = (month == 4)
gen byte d_m_5 = (month == 5)
gen byte d_m_6 = (month == 6)
gen byte d_m_7 = (month == 7)
gen byte d_m_8 = (month == 8)
gen byte d_m_9 = (month == 9)
gen byte d_m_10 = (month == 10)
gen byte d_m_11 = (month == 11)
gen byte d_m_12 = (month == 12)

gen d_y_2008 = (year == 2008)
gen d_y_2009 = (year == 2009)
gen d_y_2010 = (year == 2010)
gen d_y_2011 = (year == 2011)
gen d_y_2012 = (year == 2012)
gen d_y_2013 = (year == 2013)
gen d_y_2014 = (year == 2014)

// generating income  x liquor group interactions
foreach inc of varlist inc_* {
	foreach group of varlist d_g_*{
		gen byte `inc'_x_`group' = `inc'*`group'
	}
}


// generating proof x income group
gen proof_x_inc_lt30 = proof*inc_lt30
gen proof_x_inc_30_to_70 = proof*inc_30_to_70
gen proof_x_inc_70_to_100 = proof*inc_70_to_100
gen proof_x_inc_100_plus = proof*inc_100_plus

// generating imported x income group
gen byte imported_x_inc_lt30 = imported*inc_lt30
gen byte imported_x_inc_30_to_70 = imported*inc_30_to_70
gen byte imported_x_inc_70_to_100 = imported*inc_70_to_100
gen byte imported_x_inc_100_plus = imported*inc_100_plus

//// sampling cases
/*
gen sort_order = _n
egen select = tag(case)
gen rand = runiform()
sort select rand
replace select = (_n > (_N-5000)) // define number of cases here
bys case (select): replace select = select[_N]
sort sort_order
drop sort_order rand
keep if select == 1
*/
/////

//// sampling within case to further reduce sample size
/*
sample 10 if choice == 0 & product != 0, by(case) count
*/


//generating product dummies
foreach i of numlist 1/101{
	gen byte d_p_`i' = (product == `i')
}


gen byte holiday = d_m_11 + d_m_12 + d_m_1 // Nov or Dec or Jan

gen byte d_q_1 = d_m_1 + d_m_2 + d_m_3
gen byte d_q_2 = d_m_4 + d_m_5 + d_m_6
gen byte d_q_3 = d_m_7 + d_m_8 + d_m_9
gen byte d_q_4 = d_m_10 + d_m_11 + d_m_12

//export delimited asclogit_test.csv, replace

replace product = 102 if product == 101 & d_g_vod == 1
replace product = 103 if product == 101 & d_g_sch == 1
replace product = 104 if product == 101 & d_g_brb == 1
replace product = 105 if product == 101 & d_g_whs == 1
replace product = 106 if product == 101 & d_g_teq == 1
replace product = 107 if product == 101 & d_g_otr == 1

clogit choice p_x_* imported_x_* proof_x_* d_g_gin d_g_vod d_g_rum d_g_sch d_p_1-d_p_101, group(case) iter(150)
/*
mixlogit choice imported_x_* proof_x_* d_g_gin d_g_vod d_g_rum d_g_sch d_p_1-d_p_101, group(case) rand(p_x_*) iter(150) nrep(1)
*/



/* Adjusting constants to account for choice based sample */

predict xb, xb

gen adj_util = xb - ln(sample_weight/pop_weight)

gen num = exp(adj_util)
egen denom = total(exp(adj_util)), by(case)
gen prob = num/denom

local alpha1 = _b[p_x_inc_lt30]
local alpha2 =  _b[p_x_inc_30_to_70]
local alpha3 = _b[p_x_inc_70_to_100]
local alpha4 = _b[p_x_inc_100_plus]

gen ed_ij = (`alpha1'*inc_lt30 + `alpha2'*inc_30_to_70 + `alpha3'*inc_70_to_100 + `alpha4'*inc_100_plus)*price*(1-prob)

egen ed_j = mean(ed_ij), by(product)

sum ed_j, det


log close
