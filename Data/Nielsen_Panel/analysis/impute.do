capture log close
clear

use "../../Data/Nielsen_Panel/analysis/blp_sample.dta"
gen date_q = qofd(date_n)
format date_q %tq

fillin mkt product // generating record for all products in all markets

/* Filling in data. Fillin leaves many missings, but they can be replaced easily */

sort mkt M

foreach var of varlist M date_q d_c_* {
	by mkt: replace `var' = `var'[1]
}

sort product mkt

foreach var of varlist d_gin d_vod d_sch d_whs d_brb d_teq d_otr d_rum d_s_* import proof{
	bys product (`var'): replace `var' = `var'[1]
}

egen med_price = median(price), by(product)

replace avg_price = med_price if avg_price == .

/* J = 193 */

replace share = 1 / (M + 1 + 193) if share == .

/* Imputation and Laplace adjustment over */

gen delta = log(share) - log(outside_share) // explanatory variable
gen cons = 1
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
	
