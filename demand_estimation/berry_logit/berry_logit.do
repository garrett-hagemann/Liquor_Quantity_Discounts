capture log close
clear
version 14
set more off
log using "berry_logit.txt", text replace

/* FLAGS */
local price_sched_data = 0

/* This file estimates the so-called Berry-logit. This is basically
a version of BLP without random coefficients. This file makes the
relevant instruments:
	exogenous variable squares, interactions of exo vars,
	"other product" characteristics, etc.
*/

use "../../Data/combined/merged_blp_sample"

drop if purchase_year == 2008 // only have partial month here

gen delta = log(share) - log(outside_share) // explanatory variable
gen cons = 1



/* Making instruments */
if `price_sched_data'{
	gen blp_inst1 = proof_alcohol_cont^2
	gen blp_inst2 = price_per_case^2
}
local i = 3
foreach var in  d_gin d_vod d_rum d_sch d_whs d_brb {
	egen blp_inst`i' = total(`var'), by(mkt)
	replace blp_inst`i' = blp_inst`i' - `var'
	local i = `i' + 1
}

/* making month and year dummies for easier export */
levelsof purchase_month, local(months)
foreach month of local months{
	gen d_m`month' = (purchase_month == `month')
}

gen holiday = (d_m11 + d_m12)

levelsof purchase_year, local(years)
foreach year of local years{
	gen d_y`year' = (purchase_year == `year')
}

if `price_sched_data'{
	ivregress gmm delta d_gin d_vod d_rum d_sch d_brb d_whs proof ///
		 (price = price_per_case ), nocons first
 }
else{	 
	ivregress gmm delta d_gin d_vod d_rum d_sch d_brb d_whs holiday ///
		(price = blp_inst* d_m* d_y* cons ), nocons first 
}

keep if e(sample)

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


/* Calculating the "retail type" associated with each (p,q) bundle. This treats
each (p,q) as if it came from a monopolist facing independent demand and independent cost.
This is clearly at odds with the model that was estimated that was just estimated. However,
to do it correctly, each (p,q) must reliably be assocated with the same retailer. This is
difficult due to the Nielsen data not identifying retailers beyond a catch-all flag for
many of the products. However, this can serve as a first step. In particular, the logit demand
assumes everything is a substitute so we can bound the misspecification perhaps. */

local alpha = _b[price]
gen markup = (1/`alpha')*(1 + exp(delta))
/*gen rho = 
gen lambda = price - markup - rho */


predict xi, resid

save berry_logit, replace
export delimited using "berry_logit.csv", replace delim(",") nolab

log close
