// Using Fake Data
clear
set seed 1987
set obs 100000

// iding consumers
gen id = _n

// generating markets
gen mkt = mod(_n, 4)

// disturbing prices a little
gen p1 = .
replace p1 = 17 + 0*rnormal(0,.1) if mkt == 0
replace p1 = 22 + 0*rnormal(0,.1) if mkt == 1
replace p1 = 14 + 0*rnormal(0,.1) if mkt == 2
replace p1 = 20 + 0*rnormal(0,.1) if mkt == 3

gen p2 = .
replace p2 = 5 + 0*rnormal(0.1) if mkt == 0
replace p2 = 10 + 0*rnormal(0.1) if mkt == 1
replace p2 = 12 + 0*rnormal(0.1) if mkt == 2
replace p2 = 8 + 0*rnormal(0.1) if mkt == 3

// utility params
local b = -.2
local c = 0

// logit shares
gen s1 = exp(`b'*p1 + `c') / (1 + exp(`b'*p1 + `c') + exp(`b'*p2 + `c'))
gen s2 = exp(`b'*p2 + `c') / (1 + exp(`b'*p1 + `c') + exp(`b'*p2 + `c'))
gen s0 = 1 - s1 - s2

gen t1 = s1
gen t2 = s1+s2

// log shares
gen log_s1 = ln(s1)
gen log_s2 = ln(s2)
gen log_s0 = ln(s0)

// subtracting off outside optoin
gen y1 = log_s1 - log_s0
gen y2 = log_s2 - log_s0

// generating binary choice for logit
gen rand = runiform()
gen choice = 0
replace choice = 1 if rand <= t1
replace choice = 2 if rand <= t2 & rand > t1

gen p = 0
replace p = p1 if choice == 1
replace p = p2 if choice == 2

/* forming into form for asclogit */
gen p0 = 0

rename p p_paid

reshape long p s t log_s y, i(mkt rand choice p_paid id) j(option)
gen chosen = (choice == option)

asclogit chosen p, case(id) alt(option) nocons // need to know what each alt is for each id

// regressions
keep if chosen == 1 // getting rid of extra obs so code below works as intended
egen choice_totals = count(mkt), by(mkt choice)
egen total = count(mkt), by(mkt)

gen s_hat = choice_totals / total
gen log_s_hat = ln(s_hat)
egen inside_share = total((choice != 0))
replace inside_share = inside_share / 100000

gen norm_log_s_hat = log_s_hat - ln(1-inside_share)

//preserve
	keep if choice != 0
	reg norm_log_s_hat p, nocons cluster(mkt) // may need some sort of weighting for std err
//restore

end
// using some real data
clear
webuse lbw
logit low age lwt i.race smoke ptl ht ui

/* all one market here. so just need to get s_hat. But linear approach
won't work unless there is variation in "market share." So doesn't work
here unless you observe multiple time periods / geographic locations */

egen choice_total = total(low)
egen total = count(id)

gen s_hat = choice_total / total

reg s_hat age lwt i.race smoke ptl ht ui
