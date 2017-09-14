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

gen delta_w_profit = lin_w_profit - base_w_profit
gen delta_r_profit = lin_r_profit - base_r_profit
gen delta_cs = lin_cs - base_cs

gen total_change = delta_w_profit + delta_r_profit + delta_cs
gen pos_change = (total_change > 0) if total_change != .
gen q_weight = exp(-Q_mode)
gen p_diff = base_wavg_p - lin_avg_p
egen total_options = rownonmiss(actual_p*) if _merge_ps == 3 & !no_steps
gen big_change = (abs(total_change) > 10000) if total_change != .
gen good = (_merge_ps == 3) & !big_change & !pos_change

tab total_options
tab pos_change

// descriptive stats for estimted parameters

tabstat c_mode b_mode zeta_mid, statistics(mean sd p25 p50 p75) columns(statistics)
tabstat delta_w_profit delta_r_profit delta_cs total_change, statistics(mean sd p5 p25 p50 p75 p95) columns(statistics)

kdensity c_mode, title("") ///
	xtitle("Wholesaler Marginal Cost")
graph export c_density.pdf, replace

kdensity b_mode, bwidth(2) title("") ///
	xtitle("b in Beta(1,b)")
graph export b_density.pdf, replace 

kdensity b_mode if b_mode < 39, bwidth(2)title("") ///
	xtitle("b in Beta(1,b)")
graph export b_mode_less.pdf, replace 

kdensity zeta_mid if zeta_mid <= 50, bwidth(1.0) ///
	title("") ///
	xtitle("Cost of Additional Option (Midpoint)")
graph export zeta_density.pdf, replace

twoway (scatter lin_avg_p base_wavg_p if pos_change) ///
	(scatter lin_avg_p base_wavg_p if !pos_change) ///
	(function y=x, range(10 60)), ///
	xtitle("With Quantity Discounts") ///
	ytitle("Without Quantity Discoutns") ///
	title("Average Retailer Price") ///
	legend(order(1 2) ///
		label(1 "Positive Total Change") ///
		label(2 "Negative Total Change"))
graph export r_price_comp.pdf, replace

twoway (scatter lin_avg_p base_wavg_p if good) ///
	(function y=x, range(10 60)), ///
	xtitle("With Quantity Discounts") ///
	ytitle("Without Quantity Discoutns") ///
	title("Average Retailer Price") ///
	legend(off)
graph export r_price_comp_good.pdf, replace

tabstat c_mode b_mode zeta_mid, statistics(mean sd p25 p50 p75) columns(statistics) by(pos_change)
tabstat delta_w_profit delta_r_profit delta_cs total_change, statistics(mean sd p5 p25 p50 p75 p95) columns(statistics) by(pos_change)

graph box total_change if good, over(total_options) noout ///
	ytitle("Total Welfare Change (Monthly $)") ///
	title("Total Welfare Change by Number of Options")
graph export change_by_parts.pdf,replace

twoway (scatter total_change c_mode) ///
	(lowess total_change c_mode, bwidth(.5) lwidth(thick)) if good, ///
	ytitle("Total Welfare Change") xtitle("Wholesaler Marginal Cost") //
graph export welf_vs_c.pdf, replace

twoway (lowess total_change c_mode, bwidth(.5) lwidth(thick)) if good, ///
ytitle("Total Welfare Change") xtitle("Wholesaler Marginal Cost") yline(0) //
graph export welf_vs_c_no_scatter.pdf, replace

twoway (scatter total_change b_mode) ///
	(lowess total_change b_mode, bwidth(.5) lwidth(thick)) if good, ///
	ytitle("Total Welfare Change") xtitle("b in Beta(1,b)")
graph export welf_vs_b.pdf, replace

twoway (lowess total_change b_mode, bwidth(.5) lwidth(thick)) if good, ///
ytitle("Total Welfare Change") xtitle("b in Beta(1,b)") yline(0) //
graph export welf_vs_b_no_scatter.pdf, replace

twoway (lowess delta_w_profit b_mode) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("b in Beta(1,b)")
graph export w_profit_vs_b_no_scatter.pdf, replace

twoway (lowess delta_r_profit b_mode) if good, ///
	ytitle("Change in Average Retailer Profits") ///
	xtitle("b in Beta(1,b)")
graph export r_profit_vs_b_no_scatter.pdf, replace

twoway (lowess delta_w_profit c_mode) if good, ///
	ytitle("Change in Wholesaler Profits") ///
	xtitle("Wholesaler Marginal Cost")
graph export w_profit_vs_c_no_scatter.pdf, replace

twoway (lowess delta_r_profit c_mode) if good, ///
	ytitle("Change in Average Retailer Profits") ///
	xtitle("Wholesaler Marginal Cost")
graph export r_profit_vs_c_no_scatter.pdf, replace

twoway (lowess delta_cs b_mode) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("b in Beta(1,b)")
graph export cs_vs_b_no_scatter.pdf, replace

twoway (lowess delta_cs c_mode) if good, ///
	ytitle("Change in Consumer Surplus") ///
	xtitle("Wholesaler Marginal Cost")
graph export cs_vs_c_no_scatter.pdf, replace

graph box delta_w_profit if good, over(total_options) noout ///
	ytitle("Change in Wholesaler Profits") ///
	title("Wholesaler Profit Change by Number of Options")
graph export w_change_by_parts.pdf,replace

graph box delta_r_profit if good, over(total_options) noout ///
	ytitle("Change in Average Retailer Profits") ///
	title("Change in Average Retailer Profits by Number of Options")
graph export r_change_by_parts.pdf,replace

graph box delta_cs if good, over(total_options) noout ///
	ytitle("Change in Consumer Surplus") ///
	title("Change in Consumer Surplus by Number of Options")
graph export cs_change_by_parts.pdf,replace

/* becuase I forgot to flip around the retailer profit calculations...*/
foreach var of varlist r_t_prof*{
	replace `var' = -`var'
}

preserve
	keep if product == 64 & purchase_year == 2011 & purchase_month == 5 // Bacardi Rum
	reshape long t_cut_ r_t_prof_, i(product purchase_month purchase_year) j(option) string
	//replace t_cut_ = 1.0 if option == "top"
	twoway (connected r_t_prof_ t_cut) if option != "top", ///
		title("Bacardi Superior 750ML, May 2011") ///
		ytitle("Change in Type {&lambda} Retailer Profit") ///
		xtitle("Type {&lambda}") ///
		note("Approximately 38% of the distribution in first segment and 34% in second segment")
	graph export r_change_by_type.pdf, replace
restore

log close
