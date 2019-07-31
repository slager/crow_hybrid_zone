library(magrittr)
library(plyr)
library(dplyr)

read.csv("parental define from Clumpp 0.98.csv",stringsAsFactors=F) -> structure
read.csv("TableS1_crow_labels_pub_molecolSEX.csv",stringsAsFactors=F) -> sex
names(sex)[which(names(sex)=="ddrad.id")] <- "id"

sex$id[which(sex$id=="nca01")] <- "ca01"
sex$id[which(sex$id=="nca02")] <- "ca02"
sex$id[which(sex$id=="nca03")] <- "ca03"
sex$id[which(sex$id=="nca04")] <- "ca04"

merge(x=structure,y=sex,by="id") -> df

t.test(formula = A ~ sex, data=df,var.equal=F)
table(df$sex) # Number of males and females in sample (4 NAs)
