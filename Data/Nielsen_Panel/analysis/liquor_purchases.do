capture log close
clear
set more off
log using "liquor_purchases.txt", text replace

display "$S_TIME  $S_DATE"

/* Notes:

1) Files are huge. Importing is really slow. Patience, friend.

2) Make sure the codes are set at the begining. These are an attempt to speed
up execution. This avoids rebuilding data files if need be.
*/

local rebuild_upcs = 1 //rebuild UPC file from scratch
local rebuild_panelists = 1 // rebuild panelist and their associated trips

local years = "2008 2009 2010 2011 2012 2013 2014"

if `rebuild_upcs' {
	/* Importing UPCs to keep just the liquor ones */
	import delimited using "../nielsen_extracts/HMS/Master_Files/Latest/products.tsv", stringcol(1)

	/* Keeping only Liquor records. CAn narrow it down to Alcoholic beverage first,
	but no need. Only other Alocohol is coded as Beer or Wine. There may be some
	marginal cases that differ between the NYS data and the Nielsen Data. Right now, 
	just working with the Liquor records for simplicity */
	keep if product_group_desc == "LIQUOR"  // corresponding product_group_code is 5002

	/* Uncomment the following line to remove multi-packs. This might make matching
	brand strings easier between the NYS and Nielsen datasets. The price schedules do
	contain multipacks, but it seems more complicated to model. The multipacks represent
	a form of retail quantity discounts. The focus here is on wholesale quantity discoutns

	drop if multi >= 1 */

	save liquor_upcs, replace


	clear
}

foreach year in `years' {
	if `rebuild_panelists' {
		/* Import Panelists and keep only NY records */
		
		import delimited using "../nielsen_extracts/HMS/`year'/Annual_Files/panelists_`year'.tsv", clear
		
		keep if fips_state_desc == "NY" // associated code is 36
		
		save ny_panelists_`year', replace
		
		clear
		
		import delimited using "../nielsen_extracts/HMS/`year'/Annual_Files/trips_`year'.tsv"
		
		merge m:1 household_code panel_year using ny_panelists_`year', gen(_merge_panelists) keep(3)
		
		save ny_trips_`year', replace
	}

	/* Importing purchases. There are many purchases to one trip. For now, we'll
	drop all other purchases at that trip. That information may be important for 
	instrumenting purposes or some sort of multi-purchase demand model */
	import delimited using "../nielsen_extracts/HMS/`year'/Annual_Files/purchases_`year'.tsv", clear

	/*Now we need to sort the purchases to get only Liquor purchases in NYS. Merging
	is slow. It is faster if the merging file is smaller. So, we're going to get just
	Liquor purchases than merge down to just the NYS records. */

	// Identifying liquor purchases
	merge m:1 upc upc_ver_uc using liquor_upcs, gen(_merge_upcs) keep(3) // m:1 because multiple purchases of same item

	// WHich liquor purchases were from trips in NY
	merge m:1 trip_code_uc using ny_trips_`year', gen(_merge_trips) keep(3) // m:1 because many purcahses to 1 trip

	save ny_liquor_purchases_`year', replace
	
	clear
}

foreach year in `years' {
	append using ny_liquor_purchases_`year'
}

gen date_n = date(purchase_date, "YMD")

format date_n %td

gen purchase_year = year(date_n)
gen purchase_month = month(date_n)
gen purchase_quarter = quarter(date_n)

egen brand_string = concat(brand_descr upc_descr size1_amount size1_units), punct(";")
egen upc_id = group(upc upc_ver_uc purchase_month purchase_year)

/* Thinning out data set to just NYC. Could possibly do earlier, but easier
to do here */
gen NYC = (fips_county_descr == "BRONX" | fips_county_descr == "NEW YORK" | fips_county_descr == "KINGS" | fips_county_descr == "QUEENS" | fips_county_descr == "RICHMOND")
keep if NYC == 1


save ny_liquor_purchases_all_years, replace

display "$S_TIME  $S_DATE"

log close
