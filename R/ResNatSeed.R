#' Computation of Mixture Suitability Index (MSI) and Mixture Reliability Index (MRI) for a mixture/donor grassland composition 
#' in a restoration site
#'
#' @description The function allows to compute the Mixture Suitability Index (MSI) and Mixture Reliability Index (MRI) for a mixture/donor grassland composition 
#' in a restoration site based on topographic features, i.e. elevation, slope and aspect.
#'
#' @param trainingDB default is NULL. When NULL the training dataset is derived from 4081 vegetation surveys 
#' carried out in Piedmont Region, Italy. A customized training dataset can be generated with \code{\link[ResNatSeed]{trainingDB}}
#' and then provided to "trainingDB" argument.
#' @param composition dataframe with the mixture/donor grassland composition. The dataframe 
#' consists of two columns: 1) species code abbreviated in CEP names format and 2) abundance 
#' of the species. When "trainingDB" is NULL, the species CEP names are retrievable from 
#' \code{\link[ResNatSeed]{cep.piem}} data. When "trainingDB" is not NULL, i.e. the user does not use
#' default settings, the training dataset to be indicated in the "trainingDB" argument must be firstly 
#' generated through \code{\link[ResNatSeed]{trainingDB}}. In this case, species CEP names that may 
#' be indicated in the mixture/donor grassland composition in the "composition" argument are retrievable from \code{\link[ResNatSeed]{trainingDB}} 
#' @param elevation Elevation of restoration site, expressed in meters above sea level (m a.s.l.)
#' @param slope Slope of restoration site, expressed in degrees (°)
#' @param aspect Aspect of restoration site, expressed in degrees from North (°N)
#' @return A list with three outputs:
#' - **DESCRIPTIVES**: dataframe with the descriptive information related to species used in the mixture/donor grassland composition:
#' \describe{
#' \item{*cep.names*}{Species name in CEP format}
#' \item{*species*}{Full species name}
#' \item{*cases.number*}{Number of cases related to each species}
#' \item{*min.ele*}{Minimum elevation at which a species has been detected in the training dataset}
#' \item{*max.ele*}{Maximum elevation at which a species has been detected in the training dataset}
#' \item{*min.slope*}{Minimum slope at which a species has been detected in the training dataset}
#' \item{*max.slope*}{Maximum slope at which a species has been detected in the training dataset}
#' \item{*min.south*}{Minimum southness value at which a species has been detected in the training dataset}
#' \item{*max.south*}{Maximum southness value at which a species has been detected in the training dataset}  
#' }
#' - **SURVEYED AND PREDICTED ABUNDANCE**: dataframe with the following information:
#' \describe{
#' \item{*cep.names*}{Species name in CEP format}
#' \item{*species*}{Full species name}
#' \item{*predicted.Abundance*}{Abundance of each species predicted by the best Generalized Additive Model (GAM) based on elevation, slope, and aspect of restoration site provided by the user}
#' \item{*predicted.AbundanceMax*}{Maximum abundance achievable by a species under the best possible ecological conditions, identified by all possible combinations of altitude, slope and southness. This maximum abundance is predicted by the best GAM and the limits of the three topographic variables are indicated in the output of the descriptives }
#' \item{*ratio*}{Ratio between the predicted abundance and the maximum abundance achievable. This ratio indicates how far (ratio = 0) or close (ratio = 1) a species is from its ecological optimum. }
#' \item{*R2.adj*}{R square adjusted of the best GAM}
#' \item{*RMSE*}{Root Mean Square Error of the best GAM}
#' \item{*Mixture.abundance*}{Abundance of a species indicated in the mixture/donor grassland composition}
#' \item{*Expected.abundance*}{Expected abundance (the most achievable) of a species in a restoration site with the topographical features provided by the user}  
#' }
#' - **INDEXES**: \describe{
#' \item{*MSI*}{Mixture Suitability Index (MSI). The MSI ranges from 0 (bad) to 1 (optimal). This index ...}
#' \item{*MRI*}{Mixture Reliability Index (MRI). The MRI ranges from 0 (bad) to 1 (optimal). This index ...} 
#' }
#' @import vegan tidyverse mgcv performance dplyr
#' @export
#' 
#' 

ResNatSeed<-function(trainingDB=NULL,composition,elevation,slope,aspect){
  
  #if no reference database (obtained from the 'trainingB' function) is indicated, the
  #Piedmontese database will be used
  if (is.null(trainingDB)){
    source<-env.veg.sel.new
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
                                        cases.number=nrow(merge.selected.unique),
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
                                     predicted.Abundance=species.predicted,
                                     predicted.AbundanceMax=max(species.predicted.max))
    
    df.species.predicted$ratio<-round(df.species.predicted$predicted.Abundance/df.species.predicted$predicted.AbundanceMax,2)
    
    
    if (exists("def")==FALSE){
      df.species.predicted$R2.adj<-NA
      df.species.predicted$RMSE<-NA
    }else {
      
      if (is.na(df.species.predicted$predicted.AbundanceMax)){
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
  colnames(db)[8]<-"Mixture.abundance"
  colnames(db)[1]<-"cep.names"
  db$Expected.abundance<-round(db$Mixture.abundance*db$ratio,2)
  db1<-db[complete.cases(db),]#delete NA
  
  average<-list()
  average<-data.frame(
    MSI=round(sum(db1$Expected.abundance)/sum(db1$Mixture.abundance),2),#Mixture.Suitability.index
    MRI=round(sum(db1$Mixture.abundance)/sum(composition$sra),2))#Model.reliability.index
  
  table<-list()
  table<-db
  
  table.descriptive<-do.call(rbind,descriptive)#merging the tables of the predicted values
  
  return(list(DESCRIPTIVES=table.descriptive,
              SURVEYED_AND_PREDICTED_ABUNDANCE=table,
              INDEXES=average))
  
}
