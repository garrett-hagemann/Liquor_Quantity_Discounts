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
gen byte inc_70plus = (household_income >= 26)
//gen byte inc_100_plus = (household_income > 26)

gen inc = .
replace inc = 2500 if household_income == 3
replace inc = 6500 if household_income == 4
replace inc = 9000 if household_income == 6
replace inc = 11000 if household_income == 8
replace inc = 13000 if household_income == 10
replace inc = 17500 if household_income == 11
replace inc = 22500 if household_income == 13
replace inc = 27500 if household_income == 15
replace inc = 32500 if household_income == 16
replace inc = 37500 if household_income == 17
replace inc = 42500 if household_income == 18
replace inc = 47500 if household_income == 19
replace inc = 55000 if household_income == 21
replace inc = 65000 if household_income == 23
replace inc = 85000 if household_income == 26
replace inc = 150000 if household_income == 27


gen p_x_inc_lt30 = price*inc_lt30
gen p_x_inc_30_to_70 = price*inc_30_to_70
gen p_x_inc_70plus = price*inc_70plus

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
gen proof_x_inc_70plus = proof*inc_70plus
gen proof_x_white = proof*white


// generating imported x income group
gen byte imported_x_inc_lt30 = imported*inc_lt30
gen byte imported_x_inc_30_to_70 = imported*inc_30_to_70
gen byte imported_x_inc_70plus = imported*inc_70plus
gen byte imported_x_white = imported*white

// white x size
gen byte w_x_750ML = white*d_s_750ML
gen byte w_x_1L = white*d_s_1L
gen byte w_x_175L = white*d_s_175L

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
foreach i of numlist 1/108{
	gen byte d_p_`i' = (product == `i')
}


gen byte holiday = d_m_11 + d_m_12 + d_m_1 // Nov or Dec or Jan

gen byte d_q_1 = d_m_1 + d_m_2 + d_m_3
gen byte d_q_2 = d_m_4 + d_m_5 + d_m_6
gen byte d_q_3 = d_m_7 + d_m_8 + d_m_9
gen byte d_q_4 = d_m_10 + d_m_11 + d_m_12

//export delimited asclogit_test.csv, replace


gen lip = ln(inc - price)
gen lip_x_inc_lt30 = lip*inc_lt30
gen lip_x_inc_30_to_70 = lip*inc_30_to_70
gen lip_x_inc_70plus = lip*inc_70plus

nlogitgen type = product(oo: 0, gin: 27 | 41 | 101, vod: 1 | 3 | 6 | 8 | 11 | 13 | 14 | 15 | 17 | 19 | 20 | ///
	22 | 25 | 31 | 38 | 45 | 46 | 47 | 49 | 51 | 52 | 54 | 57 | 59 | 60 | 66 | 79 | ///
	82 | 83 | 85 | 87 | 93 | 96 | 97 | 102, rum: 2 | 16 | 18 | 21 | 24 | 28 | 29 | 34 | ///
	37 | 48 | 50 | 53 | 55 | 58 | 67 | 68 | 70 | 71 | 76 | 80 | 84 | 89 | 100 | 108, ///
	sch: 30 | 33 | 61 | 90 | 92 | 94 | 103, brb: 7 | 23 | 39 | 42| 73 | 104, ///
	whs: 5 | 56 | 91 | 105, teq: 32 | 35 | 65 | 106, otr: 4 | 9 | 10 | 12 | 26 | 36 | ///
	40 | 43 | 44 | 62 | 63 | 64 | 69 | 72 | 74 | 75 | 77 | 78 | 81 | 86 | 88 | 95 | 98 | 99 | 107)


keep if year == 2011
clogit choice lip_x_inc_* imported_x_* proof_x_* d_s_* d_p_*, group(case) iter(150)
/*
clogit choice lip_x_* imported_x_* proof_x_* d_s_750ML d_s_1L d_s_175L w_x_* d_p_*, group(case) iter(150)
asclogit choice lip_x_* imported_x_* proof_x_* d_s_750ML d_s_1L d_s_175L w_x_*, case(case) alt(product)
matrix tmp_b0 = e(b)
matrix tmp_tau = (1,1,1,1,1,1,1,1,1)
matrix rownames tmp_tau = y1
matrix colnames tmp_tau = oo_tau:_cons gin_tau:_cons vod_tau:_cons rum_tau:_cons ///
	sch_tau:_cons brb_tau:_cons whs_tau:_cons teq_tau:_cons otr_tau:_cons
matrix b0 = tmp_b0,tmp_tau

nlogit choice lip_x_* imported_x_* proof_x_* d_s_750ML d_s_1L d_s_175L w_x_* || ///
	type: || product:, case(case) iter(200) from(b0, skip)
*/
/*
asclogit chosen cost distance rating, case(family_id) alt(restaurant) iter(100)
matrix tmp_b0 = e(b)
matrix tmp_tau = (1, 1, 1)
matrix rownames tmp_tau = y1
matrix colnames tmp_tau = fast_tau:_cons family_tau:_cons fancy_tau:_cons
matrix list tmp_tau
matrix b0 = tmp_b0,tmp_tau

nlogit chosen cost distance rating || type:  || restaurant:, case(family_id) iter(50) from(b0)
*/

//nlogit choice p_x_* imported_x_* proof_x_* || type: cons || product: , case(case) iter(200)
/*
mixlogit choice imported_x_* proof_x_* d_g_gin d_g_vod d_g_rum d_g_sch d_p_1-d_p_101, group(case) rand(p_x_*) iter(150) nrep(1)
*/

/*

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
*/

log close
