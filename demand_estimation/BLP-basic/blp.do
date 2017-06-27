capture log close
clear
set more off
log using "blp.txt", text replace

use "../berry_logit/berry_logit.dta"

blp share proof_alcohol_cont d_*, ///
	endog(price = price_per_case blp_inst1 blp_inst2 blp_inst3 blp_inst4 blp_inst5 blp_inst6 blp_inst7 blp_inst8) ///
	stochastic(price) markets(mkt)

log close
