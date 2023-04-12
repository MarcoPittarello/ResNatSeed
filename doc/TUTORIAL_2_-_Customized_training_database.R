## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
#system.file("extdata", package = "ResNatSeed")
load("~/Pacchetti R/ResNatSeed/inst/extdata/veg.composition.example.RData")
library(tidyr)

## ----setup--------------------------------------------------------------------
library(ResNatSeed)

## ----echo=FALSE---------------------------------------------------------------
veg.composition.example

## -----------------------------------------------------------------------------
training.custom<-trainingDB(data=veg.composition.example,
           spe.freq = 30,# only species occurring at least in 30 surveys will be retained
           min.spe.abundance = 1# only species with at least 1% of relative abundance will be retained
           )


## -----------------------------------------------------------------------------
training.custom$cep.names

## -----------------------------------------------------------------------------
head(training.custom$trainingDB.ResNatSeed)

## -----------------------------------------------------------------------------
donor.composition<-data.frame(
  species=c("Dactglom","Festaggr","Thymaggr","Lotucorn"),
  abundance=c(25,35,15,8)
)

donor.composition

## ----echo=TRUE----------------------------------------------------------------
RestInd(trainingDB = training.custom$trainingDB.ResNatSeed, 
        composition=donor.composition,  
        elevation=1600, 
        slope=10, 
        aspect=110
        )

