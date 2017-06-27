capture log close
clear
set more off
log using "nested_logit.txt", text replace

/* This file estimates the a nested logit with aggregate data. The nests are
defined by liquor type dummies. All necessary instruments are created within
this file.
*/

use "../../Data/Nielsen_Panel/analysis/blp_sample.dta"

drop if purchase_year == 2007 // only partial data and no price schedules

gen delta = log(share) - log(outside_share) // explanatory variable
gen cons = 1
gen date_q = qofd(date_n)
format date_q %tq


// generating group shares. Liquor type dummies can serve as group identifiers
egen tmp_gin_gs = total(share) if d_gin == 1, by(mkt)
egen tmp_vod_gs = total(share) if d_vod == 1, by(mkt)
egen tmp_rum_gs = total(share) if d_rum == 1, by(mkt)
egen tmp_sch_gs = total(share) if d_sch == 1, by(mkt)
egen tmp_brb_gs = total(share) if d_brb == 1, by(mkt)
egen tmp_whs_gs = total(share) if d_whs == 1, by(mkt)
egen tmp_teq_gs = total(share) if d_teq == 1, by(mkt)
egen tmp_otr_gs = total(share) if d_otr == 1, by(mkt)

egen group_share = rowtotal(tmp_*_gs)

gen group_delta = log(share) - log(group_share)


/* Making instruments */
// Number of other products in group in market
foreach var in gin vod rum sch whs brb teq otr{
	egen other_`var' = total(d_`var'), by(mkt)
	replace other_`var' = (other_`var' - d_`var')
}

// total number of products in market
egen total_products = count(product), by(mkt)

// average price of product across markets
egen otr_mkt_avg_p = total(avg_price), by(product)
egen total_mkts_prod = count(mkt), by(product)
replace otr_mkt_avg_p = (otr_mkt_avg_p - avg_p)/ (total_mkts_prod - 1)



/* making month and year dummies for easier export */
levelsof purchase_quarter, local(quarters)
foreach quarter of local quarters{
	gen d_q`quarter' = (purchase_quarter == `quarter')
}

levelsof purchase_year, local(years)
foreach year of local years{
	gen d_y`year' = (purchase_year == `year')
}

/* Nested Logit */
ivregress gmm delta cons d_gin d_vod d_rum d_sch d_brb d_whs d_teq imported proof ///
	d_s_750 d_s_1L d_s_175L d_c_bronx d_c_kings d_c_ny d_c_queens i.date_q ///
	(avg_price = other_* otr_mkt_avg_p), nocons
	
/* Simple Nested Logit */
ivregress gmm delta cons d_gin d_vod d_rum d_sch d_brb d_whs d_teq imported proof ///
	d_s_750 d_s_1L d_s_175L d_c_bronx d_c_kings d_c_ny d_c_queens i.date_q ///
	(avg_price group_delta = other_* otr_mkt_avg_p), nocons first

/* generating higher order intruments */
gen otr_mkt_avg_p2 = otr_mkt_avg_p^2
foreach var of varlist other_* { // squared other products
	gen `var'2 = `var'^2
}

/* generating interactions of instruments */
gen i1 = otr_mkt_avg_p*other_gin
gen i2 = otr_mkt_avg_p*other_vod
gen i3 = otr_mkt_avg_p*other_rum
gen i4 = otr_mkt_avg_p*other_sch
gen i5 = otr_mkt_avg_p*other_whs
gen i6 = otr_mkt_avg_p*other_brb
gen i7 = otr_mkt_avg_p*other_teq
gen i8 = otr_mkt_avg_p*other_otr
gen i9 = otr_mkt_avg_p*other_gin2
gen i10 = otr_mkt_avg_p*other_vod2
gen i11 = otr_mkt_avg_p*other_rum2
gen i12 = otr_mkt_avg_p*other_sch2
gen i13 = otr_mkt_avg_p*other_whs2
gen i14 = otr_mkt_avg_p*other_brb2
gen i15 = otr_mkt_avg_p*other_teq2
gen i16 = otr_mkt_avg_p*other_otr2
gen i17 = otr_mkt_avg_p2*other_gin
gen i18 = otr_mkt_avg_p2*other_vod
gen i19 = otr_mkt_avg_p2*other_rum
gen i20 = otr_mkt_avg_p2*other_sch
gen i21 = otr_mkt_avg_p2*other_whs
gen i22 = otr_mkt_avg_p2*other_brb
gen i23 = otr_mkt_avg_p2*other_teq
gen i24 = otr_mkt_avg_p2*other_otr
gen i25 = otr_mkt_avg_p2*other_gin2
gen i26 = otr_mkt_avg_p2*other_vod2
gen i27 = otr_mkt_avg_p2*other_rum2
gen i28 = otr_mkt_avg_p2*other_sch2
gen i29 = otr_mkt_avg_p2*other_whs2
gen i30 = otr_mkt_avg_p2*other_brb2
gen i31 = otr_mkt_avg_p2*other_teq2
gen i32 = otr_mkt_avg_p2*other_otr2
gen i33 = other_gin*other_vod
gen i34 = other_gin*other_rum
gen i35 = other_gin*other_sch
gen i36 = other_gin*other_whs
gen i37 = other_gin*other_brb
gen i38 = other_gin*other_teq
gen i39 = other_gin*other_otr
gen i40 = other_vod*other_rum
gen i41 = other_vod*other_sch
gen i42 = other_vod*other_whs
gen i43 = other_vod*other_brb
gen i44 = other_vod*other_teq
gen i45 = other_vod*other_otr
gen i46 = other_rum*other_sch
gen i47 = other_rum*other_whs
gen i48 = other_rum*other_brb
gen i49 = other_rum*other_teq
gen i50 = other_rum*other_otr
gen i51 = other_sch*other_whs
gen i52 = other_sch*other_brb
gen i53 = other_sch*other_teq
gen i54 = other_sch*other_otr
gen i55 = other_whs*other_brb
gen i56 = other_whs*other_teq
gen i57 = other_whs*other_otr
gen i58 = other_brb*other_teq
gen i59 = other_brb*other_otr
gen i60 = other_teq*other_otr
	
ivregress gmm delta cons d_gin d_vod d_rum d_sch d_brb d_whs d_teq imported proof ///
	d_s_750 d_s_1L d_s_175L d_c_bronx d_c_kings d_c_ny d_c_queens i.date_q ///
	(avg_price group_delta = otr_mkt_avg_p otr_mkt_avg_p2 i1-i60), nocons first

/*
ivregress gmm delta cons d_gin d_vod d_rum d_sch d_brb d_whs d_teq imported proof ///
	d_s_750 d_s_1L d_s_175L d_c_bronx d_c_kings d_c_ny d_c_queens i.date_q i.product ///
	(avg_price group_delta = other_* otr_mkt_avg_p), nocons first

 Model with dummies is more general, but leads to positive price coefficient? 
Removing the constant and other dummies appears to fix this*/

	
keep if e(sample) // getting rid of observations not used in estimation

/* Calculating own price elasticities for products. Using Nested specification */
local sigma = _b[group_delta]
local alpha = _b[avg_price]
gen elasticity_own = `alpha'*avg_price/(1-`sigma')*(1-`sigma'*share/group_share - (1-`sigma')*share)

hist elasticity_own, width(.1)
graph export elasticity_hist.pdf, replace


/*

/* outputting coefficients to matrix */
matrix coefs = e(b)
preserve
	clear
	svmat coefs
	rename coefs1 price
	rename coefs2 d_gin
	rename coefs3 d_vod
	rename coefs4 d_rum
	rename coefs5 d_sch
	rename coefs6 d_brb
	rename coefs7 d_whs
	rename coefs8 holiday
	/*rename coefs8 proof_alcohol_cont
	rename coefs9 d_m1
	rename coefs10 d_m2
	rename coefs11 d_m3
	rename coefs12 d_m4
	rename coefs13 d_m5
	rename coefs14 d_m6
	rename coefs15 d_m7
	rename coefs16 d_m8
	rename coefs17 d_m9
	rename coefs18 d_m10
	rename coefs19 d_m11
	rename coefs20 d_m12
	rename coefs21 d_y2009
	rename coefs22 d_y2010
	rename coefs23 d_y2011
	rename coefs24 d_y2012
	rename coefs25 d_y2013
	rename coefs26 d_y2014 */
	export delimited using "blogit_coefs.csv", replace 
restore


predict xi, resid

save berry_logit, replace
export delimited using "berry_logit.csv", replace delim(",") nolab

log close
/*
