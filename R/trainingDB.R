#' Custom ResNatSeed training database and list of plant species that can be indicated in the mixture/donor grassland composition 
#'
#' @description This function is to be used in case the user does not want to use the training database set by default in the function \code{\link[ResNatSeed]{RestInd}}, which includes vegetation surveys carried out in the Region Piedmont, Italy. \cr
#' Starting from a database containing topographic variables (elevation, slope and aspect) and plant species abundances, the function allows to generate:
#' - the list of plant species that can be indicated in the mixture/donor grassland composition with associated their species codes in CEP names format (see \code{\link[vegan]{make.cepnames}} for more details) 
#' - a training database in the format requested as input in the \code{\link[ResNatSeed]{RestInd}} function. 
#'
#' @param data a dataframe where rows are surveys and columns are species. Such a dataframe must have the columns in the following order: Survey ID, Elevation, Slope, Aspect, and then columns with plant species abundance. See @details for additional and more specific information.
#' @param min.spe.abundance value of abundance (greater than) beyond which species will be modelled. When NULL the parameter is set to 0.
#' @param spe.freq minimum frequency (equal or greater than) of species occurrence within the training database. It is advisable to set values greater than/equal to 30 to allow appropriate modelling.  
#' @return a list containing two dataframes:\cr
#' - **"cep.names"**: a two-column dataframe reporting the full species names and their abbreviations in the CEP format (see \code{\link[vegan]{make.cepnames}} for more details). These species are those that can be indicated (as CEP names) in the mixture/donor grassland composition (i.e. in the "composition" argument of \code{\link[ResNatSeed]{RestInd}}).
#' - **"trainingDB.ResNatSeed"**: a seven-column dataframe reporting the following information: Survey ID, elevation, slope, southness, species full names, species CEP names, and species abundance. This dataframe must be used in the "trainingDB" argument of \code{\link[ResNatSeed]{RestInd}}. 
#' @details The format of the dataframe required in the "data" argument must looks as follow (the heading names do not necessarily have to be the same as those indicated, but it is essential to maintain that order):\cr
#' 
#' SurveyID   |    Elevation   |    Slope   |    Aspect   |    Species1   |    Species2   |    Species3   |    Species…\cr
#'
#'- **SurveyID** should contain an alphanumeric coding
#'- **Elevation** expressed in meters above sea level (m a.s.l.)
#'- **Slope** expressed in degrees (°)
#'- **Aspect** expressed in degrees from North (°N) as it will be converted to "southness" to avoid circular variable issues in the models run with \code{\link[ResNatSeed]{RestInd}} (Chang et al., 2004)
#'- **Species** columns contains plant species abundance, either expressed as percentage, Species Relative Cover (SRA), Species Percentage Cover (%SC)(see more details regarding the last two terminologies in \code{\link[iPastoralist]{vegetation_abundance}} or \href{ https://raw.githubusercontent.com/MarcoPittarello/iPastoralist/main/image/Wrkflw_abundance_conversion.png }{here} \cr
#' 
#' @references * Chang, C., P. Lee, M. Bai, and T. Lin. 2004. Predicting the geographical distribution of plant communities in complex terrain–a case study in Fushian Experimental Forest, north- eastern Taiwan. Ecography. 27:577–588. doi:10.1111/j.0906- 7590.2004.03852.x
#' * Marco Pittarello (2022). iPastoralist: Management, conversion, and analyses of vegetation data derived from phytosociological and point-quadrat method surveys. R package version 0.0.0.9000. https://github.com/MarcoPittarello/iPastoralist.git
#' @import vegan tidyverse dplyr
#' @export

trainingDB<-function(data,spe.freq,min.spe.abundance=NULL){
  
  db<-as.data.frame(data)
  colnames(db)[1]<-"Codril"
  colnames(db)[2]<-"elevation"
  colnames(db)[3]<-"slope"
  colnames(db)[4]<-"southness"
  db$southness<-180-abs(180-db$southness)
  rownames(db)<-db[,1]
  
  
  db.species<-db[,c(5:ncol(db))]#species column selection
  
  #exclusion of occasional species, or in any case below a minimum threshold
  if (is.null(min.spe.abundance)){
    min.spe.abundance<-0
  } else   {
    min.spe.abundance<-min.spe.abundance
  }
  db.species[db.species<=min.spe.abundance] <- 0 #exclusion of occasional species
  
  db.species<-db.species[,colSums(db.species>0)>0] #elimination of species columns with sum = 0
  
  #calculation of frequency of occurrence by species
  count<-data.frame(apply(db.species,MARGIN = 2,function(x) length(which(x>0))))
  count$species<-rownames(count)
  colnames(count)[1]<-"Freq"
  count<-count[order(count$Freq,decreasing = T),]
  count$cep.names<-make.cepnames(count$species)
  
  # Selection of species with frequency >= ....
  count.sub<-subset(count,Freq>=spe.freq)
  count.sub$species<-as.factor(count.sub$species)#species list for models
  count.sub$cep.names<-make.cepnames(count.sub$species)
  count.sub$cep.names<-as.factor(count.sub$cep.names)
  
  # database reshape with all species in one column
  env.veg<-cbind(db[,1:4],db.species)
  env.veg.rshp <- pivot_longer(data = env.veg,cols = 5:ncol(env.veg))
  colnames(env.veg.rshp)[5]<-"species"
  colnames(env.veg.rshp)[6]<-"abundance"
  
  env.veg.rshp1<-merge(env.veg.rshp,count.sub,by="species")
  env.veg.rshp2<-env.veg.rshp1[,c(2:5,1,8,6)]
  
  #extraction from the reshape database (step 5) of the selected species (step 4)
  env.veg.sel<-env.veg.rshp2[env.veg.rshp2$cep.name %in% count.sub$cep.names,]
  env.veg.sel<-as.data.frame(env.veg.sel)
  
  species<-count.sub[,c("species","cep.names")]
  species<-species[order(species$cep.names,decreasing = F),]
  
  rownames(species)<-1:nrow(species)
  
  print("A list containing two dataframes has been created:")
  print("1 - cep.names: database with the list of species suitable for modeling and their species codes (CEP names) to be used for the mixture/donor grassland composition")       
  print("2 - trainingDB.ResNatSeed: database to use in the 'RestInd' function")       
  
  return(list(cep.names=species,
              trainingDB.ResNatSeed=env.veg.sel))
  
  
  
}
