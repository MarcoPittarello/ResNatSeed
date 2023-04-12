## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ResNatSeed)

## ----echo=TRUE----------------------------------------------------------------
data("cep.piem")

head(cep.piem)


## -----------------------------------------------------------------------------
donor.composition<-data.frame(
  species=c("Bromere","Bracrupe","Festaggr","Knauarve","Silevulg"),
  abundance=c(25,18,22,8,5)
)

donor.composition

## ----echo=TRUE----------------------------------------------------------------
RestInd(trainingDB = NULL, 
        composition=donor.composition,  
        elevation=2300, 
        slope=15, 
        aspect=135
        )

