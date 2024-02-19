library(msa)
library(tidyr)

setwd("C:/Users/Burne/Documents/GitHub/Bioinformatics/BioinformaticsR")

x<- read.csv("data.csv")

attach(mtcars)

#summary
head(mtcars)

#Rows and Columns Respectively
dim(mtcars)



mtcars[1,2]
