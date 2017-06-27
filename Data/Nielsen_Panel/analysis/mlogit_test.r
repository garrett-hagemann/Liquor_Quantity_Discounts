library(mlogit)
library(readr)

setwd("E:/Box/Box Sync/Snuffles Backup/gh8728/Projects/Liquor/Data/Nielsen_Panel/analysis/")

data <- read_csv("asclogit_test.csv")

data.ml <- mlogit.data(data, choice = "choice", alt.var = "product", chid.var = "case", shape = "long")

rm(data) # removing duplicate data frame

LC.ml <- mlogit(choice ~ p_x_inc_lt30 + p_x_inc_30_to_70 + p_x_inc_70_to_100 + p_x_inc_100_plus +
                  proof + imported, data.ml)

print(summary(LC.ml))