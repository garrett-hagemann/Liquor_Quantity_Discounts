capture log close
clear

log using cf_analysis.txt, text replace

use cf_sum

/* Doing some data clean up. Looks like the estimation didn't work on these
products as the algorithm didn't move from initial guess */
gen no_steps = (Q_mode == .)
foreach var of varlist c_* b_* zeta_* base_* base_avg_p base_wavg_p lin_avg_p {
	replace `var' = . if no_steps
}

// calculating number of steps in price schedule
egen total_options = rownonmiss(actual_p*) if _merge_ps == 3 & !no_steps

gen month = month(dofm(date_m))

/* 
CORRECTION
removing impact of complexity  cost here. Already subtracted off in CF,
but easier to handle it here */
replace lin_w_profit = lin_w_profit + 2*zeta_mid
replace base_w_profit = base_w_profit + total_options*zeta_mid

gen delta_w_profit = lin_w_profit - base_w_profit
gen delta_r_profit = lin_r_profit - base_r_profit
gen delta_cs = lin_cs - base_cs

gen pct_delta_w_profit = (lin_w_profit - base_w_profit)/base_w_profit
gen pct_delta_r_profit = (lin_r_profit - base_r_profit)/base_r_profit
gen pct_delta_cs = (lin_cs - base_cs)/base_cs
gen pct_total_change = (lin_w_profit + lin_r_profit + lin_cs - base_w_profit - base_r_profit - base_cs ) / (base_w_profit + base_r_profit + base_cs)

gen inv_b = 1/b_mode

gen total_change = delta_w_profit + delta_r_profit + delta_cs
gen pos_change = (total_change > 0) if total_change != .
gen q_weight = exp(-Q_mode)
gen p_diff = base_wavg_p - lin_avg_p
gen big_change = (abs(total_change) > 10000) if total_change != .
gen good = (_merge_ps == 3) & !big_change & !pos_change

tab total_options
tab pos_change

// descriptive stats for estimted parameters

tabstat c_mode b_mode inv_b zeta_mid if good, statistics(mean sd p25 p50 p75) columns(statistics)
tabstat delta_w_profit delta_r_profit delta_cs total_change if good, statistics(mean sd p5 p25 p50 p75 p95) columns(statistics)
tabstat pct_delta_w_profit pct_delta_r_profit pct_delta_cs pct_total_change if good, statistics(mean sd p5 p25 p50 p75 p95) columns(statistics)


kdensity c_mode, title("") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export c_density.pdf, replace

kdensity b_mode, bwidth(2) title("") ///
	xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export b_density.pdf, replace

kdensity inv_b, bwidth(.01) title("") ///
	xtitle("b in Beta(1,1/b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export b_density_inv.pdf, replace 

kdensity b_mode if b_mode < 39, bwidth(2)title("") ///
	xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export b_mode_less.pdf, replace 

kdensity zeta_mid if zeta_mid <= 50, bwidth(1.0) ///
	title("") ///
	xtitle("Cost of Additional Option (Midpoint)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export zeta_density.pdf, replace

twoway (scatter lin_avg_p base_wavg_p if pos_change) ///
	(scatter lin_avg_p base_wavg_p if !pos_change) ///
	(function y=x, range(10 60)), ///
	xtitle("With Quantity Discounts") ///
	ytitle("Without Quantity Discounts") ///
	title("Average Retailer Price") ///
	legend(order(1 2) ///
		label(1 "Positive Total Change") ///
		label(2 "Negative Total Change"))
graph export r_price_comp.pdf, replace

twoway (scatter lin_avg_p base_wavg_p if good) ///
	(function y=x, range(10 60)), ///
	xtitle("Average Retailer Price With Quantity Discounts") ///
	ytitle("Average Retailer PriceWithout Quantity Discounts") ///
	title("") ///
	legend(off) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export r_price_comp_good.pdf, replace

tabstat c_mode b_mode zeta_mid, statistics(mean sd p25 p50 p75) columns(statistics) by(pos_change)
tabstat c_mode inv_b zeta_mid, statistics(mean sd p25 p50 p75) columns(statistics) by(pos_change)

graph box total_change if good, over(total_options) noout ///
	ytitle("Total Welfare Change (Monthly $)") ///
	title("Total Welfare Change by Number of Options") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export change_by_parts.pdf,replace

twoway (scatter total_change c_mode) ///
	(lowess total_change c_mode, bwidth(.5) lwidth(thick)) if good, ///
	ytitle("Total Welfare Change") xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export welf_vs_c.pdf, replace

twoway (lowess total_change c_mode, bwidth(.5) lwidth(thick)) if good, ///
ytitle("Total Welfare Change") xtitle("Wholesaler Marginal Cost") yline(0) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export welf_vs_c_no_scatter.pdf, replace

twoway (scatter total_change b_mode) ///
	(lowess total_change b_mode, bwidth(.5) lwidth(thick)) if good, ///
	ytitle("Total Welfare Change") xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export welf_vs_b.pdf, replace

twoway (lowess total_change b_mode, bwidth(.5) lwidth(thick)) if good, ///
ytitle("Total Welfare Change") xtitle("b in Beta(1,b)") yline(0) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export welf_vs_b_no_scatter.pdf, replace

twoway (lowess delta_w_profit b_mode) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_profit_vs_b_no_scatter.pdf, replace

twoway (lowess delta_w_profit b_mode) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("b in Beta(1,b)") xscale(reverse) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_profit_vs_b_no_scatter_rev.pdf, replace

twoway (lowess delta_w_profit inv_b) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("b in Beta(1,1/b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_profit_vs_b_no_scatter_inv.pdf, replace

twoway (lowess delta_r_profit b_mode) if good, ///
	ytitle("Change in Average Retailer Profits") ///
	xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export r_profit_vs_b_no_scatter.pdf, replace

twoway (lowess delta_w_profit c_mode) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_profit_vs_c_no_scatter.pdf, replace

twoway (lowess delta_r_profit c_mode) if good, ///
	ytitle("Change in Average Retailer Profits") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export r_profit_vs_c_no_scatter.pdf, replace

twoway (lowess delta_cs b_mode) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("b in Beta(1,b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_vs_b_no_scatter.pdf, replace

twoway (lowess delta_cs b_mode) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("b in Beta(1,b)") xscale(reverse) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_vs_b_no_scatter_rev.pdf, replace

twoway (lowess delta_cs inv_b) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("b in Beta(1,1/b)") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_vs_b_no_scatter_inv.pdf, replace

twoway (lowess delta_cs c_mode) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_vs_c_no_scatter.pdf, replace

graph box delta_w_profit if good, over(total_options) noout ///
	ytitle("Change in Wholesaler Profits") b1title("Number of Price Segments") ///
	note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_change_by_parts.pdf,replace

graph box delta_r_profit if good, over(total_options) noout ///
	ytitle("Change in Average Retailer Profits") ///
	title("Change in Average Retailer Profits by Number of Options") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export r_change_by_parts.pdf,replace

graph box delta_cs if good, over(total_options) noout ///
	ytitle("Change in Consumer Surplus") ///
	b1title("Number of Price Segments") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_change_by_parts.pdf,replace

/* percent changes */

graph box pct_delta_w_profit if good, over(total_options) noout ///
	ytitle("Percent Change in Wholesaler Profits") ///
	title("Percent Wholesaler Profit Change by Number of Options") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_pct_change_by_parts.pdf,replace

twoway (lowess pct_delta_w_profit c_mode) if good, ///
	ytitle("Percent Change in Wholesaler Profits") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_pct_profit_vs_c_no_scatter.pdf, replace

twoway (lowess pct_delta_w_profit b_mode) if good, ///
	ytitle("Percent Change in Wholesaler Profits") ///
	xtitle("Retailer Heterogeneity") xscale(reverse) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export w_pct_profit_vs_b_no_scatter_rev.pdf, replace

twoway (lowess pct_delta_cs b_mode) if good, ///
	ytitle("Percent Change in Consumer Surplus") ///
	xtitle("Retailer Heterogeneity") xscale(reverse) note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_pct_vs_b_no_scatter_rev.pdf, replace

twoway (lowess pct_delta_cs c_mode) if good, ///
	ytitle("Percent Change in Consumer Surplus") ///
	xtitle("Wholesaler Marginal Cost") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_pct_vs_c_no_scatter.pdf, replace

graph box pct_delta_cs if good, over(total_options) noout ///
	ytitle("Percent Change in Consumer Surplus") ///
	title("Percent Change in Consumer Surplus by Number of Options") note("") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export cs_pct_change_by_parts.pdf,replace

/* FOREGONE PROFIT CALCS */

gen foregone_prof_pct = zeta_lb / (base_w_profit)*100 // want (prof(N+1) - prof(N))/prof(N). This is that.

tabstat zeta_lb if good, stat(count mean sd p25 p50 p75) by(total_options) col(stats)
tabstat foregone_prof_pct if good, stat(count mean sd p25 p50 p75) by(total_options) col(stats)

graph box zeta_lb if good, over(total_options) noout ///
	ytitle("Foregone Profits") ///
	title("Foregone Profits by Not Offering Additional Segment") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export foregone_prof_by_parts.pdf, replace

graph box foregone_prof_pct if good, over(total_options) ///
	ytitle("Percent of Current Profits") ///
	title("Foregone Profits by Not Offering Additional Segment") ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export foregone_prof_by_parts_pct.pdf, replace

reg zeta_lb proof d_g_gin-d_g_teq d_s_750-d_s_175L i.month i.total_options if good, vce(robust) cformat(%9.3f) pformat(%5.3f)
reg foregone_prof_pct proof d_g_gin-d_g_teq d_s_750-d_s_175L i.month i.total_options if good, vce(robust) cformat(%9.3f) pformat(%5.3f)


// renaming unnamed variables
egen r_t_prof_top = rowfirst(v*)

// Calculating heterogeneity effect graph for retailers
preserve
	keep if product == 64 & purchase_year == 2011 & purchase_month == 5 // Bacardi Rum
	reshape long t_cut_ r_t_prof_, i(product purchase_month purchase_year) j(option) string
	replace t_cut_ = 1.0 if option == "top"
	twoway (connected r_t_prof_ t_cut) if option != "top", ///
		title("Bacardi Superior 750ML, May 2011") ///
		ytitle("Change in Type {&lambda} Retailer Profit") ///
		xtitle("Type {&lambda}") ///
		note("Approximately 38% of the distribution in first segment and 34% in second segment")
	graph export r_change_by_type.pdf, replace
	
	drop if t_cut_ == .
	
	set obs 1000 // set total number of interpolation points here
	gen order = _n-1
	
	// second cut
	gen tag = (order == (_N-1)/3)
	replace order = (_N-1)/3 if order == 1
	replace order = 1 if tag == 1
	drop tag
	
	//third cut
	gen tag = (order == 2*(_N-1)/3)
	replace order = 2*(_N-1)/3 if order == 2
	replace order = 2 if tag == 1
	drop tag
	
	//last cut
	gen tag = (order == 3*(_N-1)/3)
	replace order = 3*(_N-1)/3 if order == 3
	replace order = 3 if tag == 1
	drop tag
	
	ipolate t_cut_ order, gen(ip_t_cut)
	ipolate r_t_prof_ ip_t_cut, gen(ip_r_t_prof)
	
	sort order
	
	gen scaled = ip_r_t_prof*betaden(1,18.66937,ip_t_cut )
	
	twoway (connected scaled ip_t_cut, msize(0)), title("") ///
		ytitle("Change in Type {&lambda} Retailer Profit") ///
		xtitle("Type {&lambda}") ///
		note("Bacardi Superior 750ML, May 2011. Scaled by estimated density. c=7.02, b=18.67") yline(0, lcolor(black)) ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
	graph export r_change_by_type_smooth.pdf, replace
	
	
restore

sum price proof imported d_s_* if good

reg c_mode d_g_gin-d_g_teq d_s_375ML d_s_750ML d_s_175L proof if good, robust

forvalues i = 0/7 {
	gen margin`i' = (actual_p`i' / bottles_per_case) / c_mode - 1
	cumul margin`i', gen(cum`i')
}

egen min_p = rowmin(actual_p*)
gen min_margin = min_p / bottles_per_case / c_mode - 1

tabstat margin* min_margin if good, stat(mean min p50 max) col(stat)
	
twoway (function y=betaden(1,1,x), 	lwidth(medthick)) ///
	(function y=betaden(1,6,x), 	lwidth(medthick)) ///
	(function y=betaden(1,.5,x), 	lwidth(medthick)), ///
	title("") ytitle("") xtitle("") ///
	legend(label(1 "Beta(1,1)") label(2 "Beta(1,6)") label(3 "Beta(1,0.5)") ///
	ring(0) bplacement(2) bmargin(large) col(1)) ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))
graph export density_comp.pdf, replace

local a = -0.4
local b = 1

twoway (function y=exp(`a'*x)/(`b' + exp(`a'*x))*1, range(0 10) lwidth(medthick)) ///
	(function y=exp(`a'*x)/(`b' + exp(`a'*x))*.75, range(0 10) lwidth(medthick)) ///
	(function y=exp(`a'*x)/(`b' + exp(`a'*x))*.5, range(0 10) lwidth(medthick)), ///
	title("") ytitle("Market Share") xtitle("Price") ///
	legend(label(1 "1.00s(p)") label(2 "0.75s(p)") label(3 "0.50s(p)") ///
	ring(0) bplacement(2) bmargin(large) col(1) off) ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))

/*
twoway (function y=-ln(1/x - `b')/`a', range(0 .5) lwidth(medthick)) ///
	(function y=-ln(1/x - `b')/`a'*.75, range(0 .5) lwidth(medthick)) ///
	(function y=-ln(1/x - `b')/`a'*.5, range(0 .5) lwidth(medthick)), ///
	title("") ytitle("Price") xtitle("Market Share") ///
	legend(label(1 "1.00s(p)") label(2 "0.75s(p)") label(3 "0.50s(p)") ///
	ring(0) bplacement(2) bmargin(large) col(1) off) ///
	graphregion(color(white)) ylabel(,nogrid) plotregion(margin(0))	
*/
	
graph export single_crossing.pdf, replace



log close
