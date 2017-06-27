/* This file cleans the new format of the NY liquor price schedules.
This constitutes the second half of 2014 and all of 2015. Once cleaned,
subsets of the data can be stacked with the old format with a little renaming.*/

capture log close
version 14
clear all
log using data_clean_new_format.txt, text replace

local months "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"


foreach m in `months'{
	disp "Parsing 2015`m'"
	import delimited using "2015/LR-`m'-2015.txt", encoding(UTF-16) clear
	tempfile tmp2015`m'
	save `tmp2015`m''
}

local 2014MonthsNew "May Jun Jul Aug Sep Oct Nov Dec"

foreach m in `2014MonthsNew'{
	disp "Parsing 2014`m'"
	import delimited using "2014/LR-`m'-2014-new-format.txt", encoding(UTF-16) clear
	tempfile tmp2014`m'
	save `tmp2014`m''
}

clear // starting with a clean dataset

foreach m in `months'{
	disp "Appending 2015`m'"
	append using `tmp2015`m''
}


foreach m in `2014MonthsNew'{
	disp "Appending 2014`m'"
	append using `tmp2014`m''
}

gen record_num = _n

save new_format_months, replace
