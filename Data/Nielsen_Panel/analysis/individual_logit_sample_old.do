/* The other data creation file focuses on identifying liquor purchases directly.
This is a useful dataset for using BLP style estimation where no information
on the outside option is needed.

However, contained in the Homescan Panel is information on non-liquor purchases.
It seems possible to estimate demand at the individual level. However, this requires
identifying what is the outside option. For now, we're identifying the outside
option as any trip that doesn't involve purchasing liquor. There's some art here:
another option would be any PURCHASE that isn't liquor. But that seems like it
makes the outside option too big - that is, too many people would choose the
outside option (especially considering that grocery trips often involve many
many items and liquor is not even an option there). This do-file constructs
that dataset.*/

capture log close
clear
set more off

log using "individual_logit_sample.txt", text replace

/* the strategy:

1) take all NY purchases and identify which are liquor purchases
2) for Trips that do not include a liquor purchase, keep only one purchase
3) for non-liquor (NL) trips: code as same choice (0)
4) for liquor trips: code choice as UPC
5) make relevant variables: proof, liquor type dummies, size dummies. These should
be set to 0 for the outside options. Note that this is equivalent to normalizing
the mean utility of the outisde option to 0.

steps 3-5 will be done in the demand estimation do file so that they can be coded
up exactly the way the estimation wants

*/

local years = "2008" // 2009 2010 2011 2012 2013 2014"

foreach year in `years'{

	import delimited using "../nielsen_extracts/HMS/`year'/Annual_Files/purchases_`year'.tsv", clear rowrange(64000000:65000000)

	// Identifying liquor purchases
	merge m:1 upc upc_ver_uc using liquor_upcs, gen(_merge_upcs) keep(1 3)
	
	// Keeping just one purchase for non-liquor trips
	egen max_merge = max(_merge_upcs), by(trip_code_uc) // largest merge code. Should be 3
	gen liquor_trip = (max_merge == 3) // trips where one of the UPCs is a liquor UPC
	drop if liquor_trip == 1 & _merge_upcs == 1 // dropping non-liquor records for liquor trips

	bys trip_code_uc: gen trip_order = _n // ordering purchases within trips. Doesn't matter if sort is stable. Just need one record for each non-liquor trip
	
	drop if liquor_trip == 0 & trip_order > 1 // dropping additional purchases from non-liquor trips
	/*
	// Getting just NYS purchases
	merge m:1 trip_code_uc using ny_trips_`year', gen(_merge_trips) keep(3) 
	*/
	tempfile obs`year'

	save `obs`year''

}

clear

foreach year in `years'{
	append using `obs`year''
}

save individual_logit_sample, replace

log close	

