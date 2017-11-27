capture log close
clear
set more off

//ssc install parmest

log using "demand_second_stage.txt", text replace

estimates use liq_demand.ster

clogit

parmest, fast

gen prod_dummy = regexm(parm,"d_p_.*")

keep if prod_dummy == 1
keep parm estimate

gen product = substr(parm,5,.)
destring product, replace

tempfile dummies
save `dummies'

end

use prod_chars_individual

merge m:1 product using `dummies', assert(1 3) keep(3)

reg estimate d_s_375 d_s_1L d_s_175L  imported proof d_g_gin-d_g_teq

log close
