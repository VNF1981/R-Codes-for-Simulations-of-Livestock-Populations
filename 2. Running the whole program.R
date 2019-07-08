rm(list=ls(all=TRUE))
ls()

setwd("e:/Genetics/1397 research materials/Codes for GitHub")
getwd()
source("e:/Genetics/1397 research materials/Codes for GitHub/1. The whole program.R")

Output <- Valpop (totqtl = 20 ,
               chrn = 2 ,
               lenchr = 1 ,
               totsnp = 100 ,
               Basenum = 100 ,
               G = 50 ,
               Hgen = 3 ,
               nRef = 300 ,
               h2 = 0.1,
               Valgen = 7 ,
               malesintensity = 0.02 ,
               femalesintensity = 0.2,
               Mating.Method = "Random.Mating")

attach(Output)
names(Output)

LDs
accmat
accmatcomp
Inbmat
geneticimp
Genetic.Variance

rm(list=ls(all=TRUE))
ls()


