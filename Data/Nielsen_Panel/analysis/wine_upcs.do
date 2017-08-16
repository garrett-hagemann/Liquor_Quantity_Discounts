capture log close
clear

import delimited using "../nielsen_extracts/HMS/Master_Files/Latest/products.tsv", stringcol(1)

keep if product_group_desc == "WINE"

keep upc upc_ver_uc

save wine_upcs, replace


