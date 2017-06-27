clear
capture log close

log using "retailer_graphs.txt", text replace

import delimited using retailer_types.csv

/* Figuring out how many times you observe a given product */
egen count = count(product), by(product)

// reversing so most observed are on top
sort count
gen rev = -_n
sort rev

// doing just one product for now
hist lambda if product == 638, width(1) kden kdenopts(bw(1)) normal normopts(lcolor(cranberry))
return list
graph export "prod638.pdf", as(pdf) replace

log close
