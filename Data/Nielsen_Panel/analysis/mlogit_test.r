library(mlogit)
library(readr)

setwd("E:/Box/Box Sync/Snuffles Backup/gh8728/Projects/Liquor/Data/Nielsen_Panel/analysis/")

data <- read_csv("test_sample.csv", col_type=cols(date_m="c",.default="d"))

data.ml <- mlogit.data(data, choice = "choice", alt.var = "product", chid.var = "case", shape = "long")

rm(data) # removing duplicate data frame

LC.ml <- mlogit(choice ~ avg_price_inter+ap_y + 
                  d_p_1 + d_p_2 + d_p_3 + d_p_4 + d_p_5 + d_p_6 + d_p_7 + d_p_8 + d_p_9 + d_p_10 +
                  d_p_11 + d_p_12 + d_p_13 + d_p_14 + d_p_15 + d_p_16 + d_p_17 + d_p_18 + d_p_19 + d_p_20 +
                  d_p_21 + d_p_22 + d_p_23 + d_p_24 + d_p_25 + d_p_26 + d_p_27 + d_p_28 + d_p_29 + d_p_30 +
                  d_p_31 + d_p_32 + d_p_33 + d_p_34 + d_p_35 + d_p_36 + d_p_37 + d_p_38 + d_p_39 + d_p_40 + 
                  d_p_41 + d_p_42 + d_p_43 + d_p_44 + d_p_45 + d_p_46 + d_p_47 + d_p_48 + d_p_49 | 0
                , data = data.ml, method="nr")
LC.nl <- mlogit(choice ~ avg_price_inter+ap_y + 
                  d_p_1 + d_p_2 + d_p_3 + d_p_4 + d_p_5 + d_p_6 + d_p_7 + d_p_8 + d_p_9 + d_p_10 +
                  d_p_11 + d_p_12 + d_p_13 + d_p_14 + d_p_15 + d_p_16 + d_p_17 + d_p_18 + d_p_19 + d_p_20 +
                  d_p_21 + d_p_22 + d_p_23 + d_p_24 + d_p_25 + d_p_26 + d_p_27 + d_p_28 + d_p_29 + d_p_30 +
                  d_p_31 + d_p_32 + d_p_33 + d_p_34 + d_p_35 + d_p_36 + d_p_37 + d_p_38 + d_p_39 + d_p_40 + 
                  d_p_41 + d_p_42 + d_p_43 + d_p_44 + d_p_45 + d_p_46 + d_p_47 + d_p_48 + d_p_49| 0
                , data = data.ml, 
                nests=list(oo=c('0'),
                           gin=c('16','25','27','40','49'),
                           vod=c('2','4','6','9','10','11','13','18','26','28','31','32','39','42','50'),
                           rum=c('5','36','44','47'),
                           sch=c('8','19','23','24','35','37'),
                           brb=c('3','14','17','20','21','22','29','33','34','38','48'),
                           whs=c('1','7','15','43','45'),
                           otr=c('12','30','41','46')), 
                un.nest.el=TRUE)

print(summary(LC.ml))