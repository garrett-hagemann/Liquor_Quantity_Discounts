capture log close
clear
set more off

ssc install parmest

log using "demand_second_stage.txt", text replace

estimates use liq_demand.ster

parmest, fast

gen prod_dummy = regexm(parm,"d_p_.*")

keep if prod_dummy == 1
keep parm estimate

gen product = substr(parm,5,.)
destring product, replace

tempfile dummies
save `dummies'

use prod_chars_individual

merge m:1 product using `dummies', assert(1 3) keep(3)

reg estimate d_s_* imported proof

log close
