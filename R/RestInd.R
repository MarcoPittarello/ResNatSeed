#' Computation of Suitability Index (SI) and Reliability Index (RI)
#'
#' @description The function allows to compute the Suitability Index (SI) and Reliability Index (RI) for a seed mixture or 
#' donor grassland composition to be used in a restoration site based on topographic features: elevation, slope and aspect.\cr
#' 
#' see vignettes for additional details and examples
#'
#' @param trainingDB dataframe with the species eligible for the statistical modelling, selected on the basis of their 
#' frequency and abundance in the vegetation and topographical variables database. Each species is associated with the 
#' values of topographic variables (elevation, slope and aspect) and its abundance. Default is NULL. When NULL the training dataset is derived from 4081 
#' vegetation surveys carried out in Piedmont Region, (North-Western Italy). A customized training dataset can be generated 
#' with \code{\link[ResNatSeed]{trainingDB}} and then provided to "trainingDB" argument.
#' @param composition dataframe with the seed mixture or donor grassland composition. It
#' consists of two columns: 1) species code abbreviated in CEP names format and 2) abundance 
#' of the species. When "trainingDB" is NULL, the species CEP names are retrievable from 
#' \code{\link[ResNatSeed]{cep.piem}} data. When "trainingDB" is not NULL, i.e. the user does not use
#' default settings, species CEP names are retrievable from \code{\link[ResNatSeed]{trainingDB}}, which
#' generates CEP names of species suitable for modeling from a customized training dataset. 
#' @param elevation Elevation of restoration site, expressed in meters above sea level (m a.s.l.)
#' @param slope Slope of restoration site, expressed in degrees (°)
#' @param aspect Aspect of restoration site, expressed in degrees from North (°N)
#' @return A list with three outputs:
#' - **DESCRIPTIVES**: dataframe with the descriptive information related to modeled plant species:
#' \describe{
#' \item{*cep.names*}{Species name in CEP format}
#' \item{*species*}{Full species name}
#' \item{*n.obs*}{Number of observations of each species selected for modeling}
#' \item{*min.ele*}{Minimum elevation at which a species occurred}
#' \item{*max.ele*}{Maximum elevation at which a species occurred}
#' \item{*min.slope*}{Minimum slope at which a species occurred}
#' \item{*max.slope*}{Maximum slope at which a species occurred}
#' \item{*min.south*}{Minimum southness value at which a species occurred}
#' \item{*max.south*}{Maximum southness value at which a species occurred}  
#' }
#' - **SPECIES ABUNDANCES**: dataframe with the following information:
#' \describe{
#' \item{*cep.names*}{Species name in CEP format}
#' \item{*species*}{Full species name}
#' \item{*PMA*}{Predicted Maximum Abundance. It is the maximum achievable abundance of a species in the restoration site, as predicted by the best model}
#' \item{*POA*}{Predicted Optimal Abundace. It is the maximum achievable abundance of a species in its optimal ecological condition, based on all possible combinations of elevation, slope and southness.}
#' \item{*ratio*}{Ratio between the *PMA* and *POA*. This ratio indicates how far (ratio = 0) or close (ratio = 1) a species is from its ecological optimum. }
#' \item{*R2.adj*}{R square adjusted of the best Generalized Additive Model}
#' \item{*RMSE*}{Root Mean Squared Error of the best Generalized Additive Model}
#' \item{*SmDgA*}{Seed mixture or Donor grassland Abundance. Abundance of a species listed in the seed mixture or donor grassland composition imported by the user and eligible for modeling}
#' \item{*EA*}{Expected Abundance. The highest achievable abundance of a species in a restoration site, based on how far the species is from the ecological optimum (i.e. computed from the multiplication of SmDgA by the ratio)}  
#' }
#' - **INDEXES**: \describe{
#' \item{*SI*}{Suitability Index (SI). Suitability of a seed mixture or donor grassland to restore a site with specific topographic characteristics. It ranges between 0 and 1. When SI=0 the restoration site is totally beyond the optimal ecological ranges of all species of the seed mixture or donor grassland, which is therefore not appropriate for the site restoration. Conversely, when SI=1 the restoration site has the optimal ecological conditions for all species of the seed mixture or donor grassland, which is therefore perfectly appropriate for the site restoration.}
#' \item{*RI*}{Reliability Index (RI). Index of the reliability of the Suitability Index (SI). The RI ranges between 0 and 1. When RI is close to 0 it means that few to none species contribute to the computation of the SI, whereas when RI is close to 1 the SI is computed with most to all species. Therefore, the higher is the RI, the most reliable is the SI. Not all the species of the seed mixture and donor grassland composition may modeled as i) they can be missing from the training database or ii) the values of the topographic factors of the restoration site are beyond their ecological ranges (e.g. if the elevation of the restoration site is 250 m and a species as an elevation range bounded between 1000 and 3000 m, such a species cannot be modeled).} 
#' }
#' @import vegan tidyverse mgcv performance dplyr
#' @examples #creation of seed mixture dataframe by retrieving species CEP names from ResNatSeed::cep.piem
#' mixture_composition<-data.frame(species=c("Festrubr","Chaehirs","Phleprat","Onobmont"),
#'                                sc=c(5,40,20,20))
#' #run function
#' ResNatSeed::RestInd(trainingDB = NULL,
#'                     composition = mixture_composition,
#'                     elevation = 1500,
#'                     aspect = 30,
#'                     slope = 5)
#' @export

RestInd<-function(trainingDB=NULL,composition,elevation,slope,aspect){
  
  #if no reference database (obtained from the 'trainingB' function) is indicated, the
  #Piedmontese database will be used
  if (is.null(trainingDB)){
    source<-trainingPie.resnatseed
  } else   {
    source<-trainingDB
  }
  
  colnames(composition)[1]<-"cep.names"
  colnames(composition)[2]<-"sra"
  
  
  #conversion aspect to southness
  southness.set<-180-abs(180-aspect)
  
  #extraction of species codes from input database
  selezione.specie.miscela<-source[source$cep.names %in% composition[,1],]
  selezione.specie.miscela$cep.names<-factor(selezione.specie.miscela$cep.names)
  
  #GAM model per each species with for cycle
  cicli<-levels(selezione.specie.miscela$cep.names)
  lista.predicted<-list()
  descriptive<-list()
  
  for (i in cicli) {
    
    species<-filter(source,cep.names==i & abundance>0) %>% droplevels()
    
    species$sc.beta<-((species$abundance*(length(species$abundance)-1)+0.5)/length(species$abundance))/100#transformation of abundance for beta family modeling 
    
    # creation and assignment of elevation, slope and southness classes to each survey
    species$cat.ele<-cut(species$elevation, 
                         seq(0,3000,50),
                         right=FALSE,
                         labels=seq(0,2950,50))
    
    species$cat.slope<-cut(species$slope, 
                           seq(0,60,5), 
                           right=FALSE,
                           labels=seq(0,55,5))
    
    species$cat.south<-cut(species$southness, 
                           seq(0,180,10), 
                           right=FALSE,
                           labels=seq(0,175,10))
    
    # selection of surveys with maximum abundance for each class for elevation, slope and southness respectively
    select.row.elevation <- species %>% 
      group_by(cat.ele) %>%
      filter(sc.beta == max(sc.beta)) 
    
    select.row.slope <- species %>% 
      group_by(cat.slope) %>%
      filter(sc.beta == max(sc.beta)) 
    
    select.row.southness <- species %>% 
      group_by(cat.south) %>%
      filter(sc.beta == max(sc.beta)) 
    
    # merging of survey selection with maximum abundance and elimination of duplicates
    merge.selected<-rbind(select.row.elevation,select.row.slope,select.row.southness)
    merge.selected.unique<-unique(merge.selected)
    
    descriptive[[i]]<-format(data.frame(row.names = FALSE,
                                        cep.names=unique(merge.selected.unique$cep.names),
                                        species=unique(merge.selected.unique$species),
                                        n.obs=nrow(merge.selected.unique),
                                        min.ele=min(merge.selected.unique$elevation),
                                        max.ele=max(merge.selected.unique$elevation),
                                        min.slope=min(merge.selected.unique$slope),
                                        max.slope=max(merge.selected.unique$slope),
                                        min.south=min(merge.selected.unique$southness),
                                        max.south=max(merge.selected.unique$southness)),justify="centre",digit=2)
    
    
    #transformation of class values into numbers (to be able to run models)
    merge.selected.unique$cat.ele<-as.numeric(levels(merge.selected.unique$cat.ele))[merge.selected.unique$cat.ele]
    merge.selected.unique$cat.slope<-as.numeric(levels(merge.selected.unique$cat.slope))[merge.selected.unique$cat.slope]
    merge.selected.unique$cat.south<-as.numeric(levels(merge.selected.unique$cat.south))[merge.selected.unique$cat.south]
    
    
    # 4 GAM models 
    a<-tryCatch(gam(sc.beta~s(elevation)+s(slope)+s(southness)
                    +ti(elevation,southness)
                    +ti(elevation,slope)
                    ,family=betar,
                    method = "REML",
                    data=merge.selected.unique), 
                error = function(e) e)
    
    b<-tryCatch(gam(sc.beta~s(elevation)+s(slope)+s(southness)
                    #+ti(elevation,southness)
                    +ti(elevation,slope)
                    ,family=betar,
                    method = "REML",
                    data=merge.selected.unique), 
                error = function(e) e)
    
    c<-tryCatch(gam(sc.beta~s(elevation)+s(slope)+s(southness)
                    +ti(elevation,southness)
                    #+ti(elevation,slope)
                    ,family=betar,
                    method = "REML",
                    data=merge.selected.unique), 
                error = function(e) e)
    
    d<-tryCatch(gam(sc.beta~s(elevation)+s(slope)+s(southness)
                    #+ti(elevation,southness)
                    #+ti(elevation,slope)
                    ,family=betar,
                    method = "REML",
                    data=merge.selected.unique), 
                error = function(e) e)
    
    
    lista<-list(a,b,c,d)#inserting models into a list
    lista.no.error<-lista[sapply(lista, length)>2]# inclusion of models correctly run in a list
    
    #if no model has correctly fitted, the dataframe of predicted values will not be considered (val 9999),
    #Otherwise, the model with the lowest AIC is selected and put on a list.
    if(length(lista.no.error)==0){
      val.settati<-data.frame(elevation=9999,
                              slope=9999,
                              southness=9999)
      suppressWarnings(rm(def))#removal of the def argument from the environment, otherwise the result of a model is overwritten even for species (later in the cycle) that cannot be modelled
      
    } else{
      minAIC<-tryCatch(which.min(sapply(lista.no.error, AIC)),error = function(e) e)
      
      def<-lista.no.error[[minAIC]]
      
      # database of values on which to make the prediction, i.e. values set by the user  
      val.settati<-data.frame(elevation=if (elevation< min(species$elevation) |elevation> max(species$elevation)  ) {9999}else {elevation},
                              slope=if (slope<5){#conditions for accepting values = 0 when the input value is <5
                                if(min(species$slope)<5) {slope} else {9999}
                              } else{ if (slope< min(species$slope) |slope> max(species$slope)  ) {9999}else {slope} },
                              southness=if (southness.set>157.5){#conditions for accepting values = at 0 or 180
                                if(max(species$southness)>157.5) {southness.set} else {9999}
                              } else if (southness.set<22.5){
                                if(min(species$southness)<22.5) {southness.set} else {9999}
                              } else{ if (southness.set< min(species$southness) |southness.set> max(species$southness)) {9999}else {southness.set}  })
      
    }
    
    # If there are values = 9999 (i.e. non-predictable values), no output will be generated, otherwise 
    # a prediction will be made based on the set values
    if(apply(val.settati,MARGIN = 1,function(x) max(x))==9999){
      species.predicted<-NA
      species.predicted.max<-NA
    } else{
      #predicted abundance for set values by user
      species.predicted<-round(predict(def,newdata = val.settati,type = "response")*100,1)
      
      #max abundance of the i species from all possible combinations of elevation, slope and southness
      newdata.all.combinations<-expand.grid(
        elevation=seq(round(min(species$elevation),0),round(max(species$elevation),0),50),
        slope=seq(round(min(species$slope),0), round(max(species$slope),0),5),
        southness=seq(round(min(species$southness),0),round(max(species$southness),0),10))
      species.predicted.max<-round(predict(def,newdata = newdata.all.combinations,type = "response")*100,1)
      newdata.all.combinations$predicted<-species.predicted.max
    }
    
    #table with species name and predicted abundance values
    df.species.predicted<-data.frame(species= unique(species$species),
                                     PMA=species.predicted,
                                     POA=max(species.predicted.max))
    
    df.species.predicted$ratio<-round(df.species.predicted$PMA/df.species.predicted$POA,2)
    
    
    if (exists("def")==FALSE){
      df.species.predicted$R2.adj<-NA
      df.species.predicted$RMSE<-NA
    }else {
      
      if (is.na(df.species.predicted$POA)){
        df.species.predicted$R2.adj<-NA
        df.species.predicted$RMSE<-NA
      } else{
        
        r2.adj<-r2(def)
        df.species.predicted$R2.adj<-round(r2.adj$R2,2)
        df.species.predicted$RMSE<-round(rmse(def)*100,1)
      }
    }
    
    lista.predicted[[i]]<-df.species.predicted
    
  }    
  
  lista.predicted.merged<-do.call(rbind,lista.predicted)#merging the tables of the predicted values
  
  db<-merge(lista.predicted.merged,composition,by.x = "row.names",by.y = "cep.names")
  colnames(db)[8]<-"SmDgA"
  colnames(db)[1]<-"cep.names"
  db$EA<-round(db$SmDgA*db$ratio,2)
  db1<-db[complete.cases(db),]#delete NA
  
  average<-list()
  average<-data.frame(
    SI=round(sum(db1$EA)/sum(db1$SmDgA),2),#Suitability.index
    RI=round(sum(db1$SmDgA)/sum(composition$sra),2))#Reliability.index
  
  table<-list()
  table<-db
  
  table.descriptive<-do.call(rbind,descriptive)#merging the tables of the predicted values
  
  return(list(DESCRIPTIVES=table.descriptive,
              SPECIES_ABUNDANCES=table,
              INDEXES=average))
  
}
