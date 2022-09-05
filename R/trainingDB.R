#' Customized ResNatSeed training database 
#'
#' @description This function creates a training database based on vegetation and topographical variables database provided
#' by the user, which must be specified in the *trainingDB* argument of 
#' the function \code{\link[ResNatSeed]{RestInd}}.\cr
#' 
#' see vignettes for additional details and examples
#'
#' @param data dataframe with vegetation and topographical variables with columns arranged as shown in @details and vignettes.
#' @param min.spe.abundance threshold of minimum abundance (greater than) of each species in each survey. When NULL the parameter is set to 0.
#' @param spe.freq threshold of minimum frequency of each species in the surveys. It is advisable to set values greater than or equal to 30 to allow appropriate statistical modeling.  
#' @return a list containing two dataframes:\cr
#' dataframe with the list of species suitable for modeling and their species codes (CEP names) to be used 
#' for the seed mixture or donor grassland composition.
#' - **"cep.names"**: a two-column dataframe reporting the full species names and their abbreviations in the
#'  CEP format (see \code{\link[vegan]{make.cepnames}} for more details). These species are those that need to be
#'   indicated (in CEP names format) in the seed mixture or donor grassland composition 
#'   (i.e. in the \code{composition} argument of \code{\link[ResNatSeed]{RestInd}}).
#' - **"trainingDB.ResNatSeed"**: a seven-column dataframe reporting the following information: Survey ID, elevation, slope, southness, species full names, species CEP names, and species abundance. This dataframe must be used in the \code{trainingDB} argument of \code{\link[ResNatSeed]{RestInd}} function. 
#' @details The format of the dataframe required in the \code{data} argument must strictly follow this column order:\cr
#' 
#' SurveyID   |    Elevation   |    Slope   |    Aspect   |    Species1   |    Species2   |    Species3   |    Species…\cr
#'
#'- **SurveyID** should contain an alphanumeric coding to uniquely identify each survey
#'- **Elevation** expressed in meters above sea level (m a.s.l.)
#'- **Slope** expressed in degrees (°)
#'- **Aspect** expressed in degrees from North (°N) as it will be converted to "southness" to avoid circular variable issues in the statistical models through the \code{\link[ResNatSeed]{RestInd}} function (Chang et al., 2004)
#'- **Species...** columns containing each plant species abundances. Abundances must be numbers bounded between 0 and 100, which can be either a species relative abundance or a species cover (*sensu* Pittarello et al., 2016; Verdinelli et al.,2022 or see \code{\link[iPastoralist]{vegetation_abundance}} and \href{ https://raw.githubusercontent.com/MarcoPittarello/iPastoralist/main/image/Wrkflw_abundance_conversion.png }{here} \cr
#' @references * Chang, C., P. Lee, M. Bai, and T. Lin. 2004. Predicting the geographical distribution of plant communities in complex terrain–a case study in Fushian Experimental Forest, north- eastern Taiwan. Ecography. 27:577–588. doi:10.1111/j.0906- 7590.2004.03852.x
#' * Pittarello, M., Probo, M., Lonati, M., Lombardi, G., 2016. Restoration of sub-alpine shrub-encroached grasslands through pastoral practices: effects on vegetation structure and botanical composition. Appl Veg Sci 19, 381–390. https://doi.org/10.1111/avsc.12222
#' * Verdinelli, M., Pittarello, M., Caria, M.C., Piga, G., Roggero, P.P., Marrosu, G.M., Arrizza, S., Fadda, M.L., Lombardi, G., Lonati, M., Nota, G., Sitzia, M., Bagella, S., 2022. Congruent responses of vascular plant and ant communities to pastoral land-use abandonment in mountain areas throughout different biogeographic regions. Ecol Process 11, 35. https://doi.org/10.1186/s13717-022-00379-9
#' * Pittarello M. (2022). iPastoralist: Management, conversion, and analyses of vegetation data derived from phytosociological and point-quadrat method surveys. R package version 0.0.0.9000. https://github.com/MarcoPittarello/iPastoralist.git
#' @import vegan tidyr
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
  print("1 - cep.names: database with the list of species suitable for modeling and their species codes (CEP names) to be used for the seed mixture or donor grassland composition")       
  print("2 - trainingDB.ResNatSeed: database to use in the 'RestInd' function")       
  
  return(list(cep.names=species,
              trainingDB.ResNatSeed=env.veg.sel))
  
}
