extractOutliers <- function( data, minPoints, minCoher=0.7, k=2, rejCrit=3, minJacc = 0.6 ){
  parameters <- list(minPoints=minPoints, minCoher=minCoher, k=k, rejCrit=rejCrit, minJacc=minJacc)

  allData <- prepareData(data)

  #Location based
  dbscan <- getDbscan( allData, parameters )
  parameters$eps <- dbscan$eps
  groupCluster <- as.integer(names(which.max(table(dbscan$cluster))))
  dataWithCluster <- bindClusterToData(dbscan, allData)

  #Data-driven
  pointsWithOCFlag <- applyPCA( dataWithCluster, parameters ) #points with oc flag; k no of components
  outlierAndNoiseFree <- getOCFreeData( pointsWithOCFlag ) #points to be included in final dataset and found by PCA
  noiseToBeIncluded <- getNoiseWithCoherLimit(pointsWithOCFlag, parameters)

  #Application-driven
  madOfVariablePerClusterDF <- calculateMadOfVariablePerCluster( pointsWithOCFlag, dbscan, parameters )
  additionalPoints <- processOCs( pointsWithOCFlag, madOfVariablePerClusterDF, parameters )
  pointsToKeep <- rbind(outlierAndNoiseFree, noiseToBeIncluded, additionalPoints)


  nonoutliers <- data[data$ID %in% row.names(pointsToKeep),]
  outliers <- data[!data$ID %in% row.names(pointsToKeep),]

  x <- list(nonoutliers = nonoutliers, outliers = outliers,parameters = parameters)

  return (x)
}

#TODO column + data validations
prepareData <- function( data ){
  inData <- data[-1]
  row.names(inData) <- data$ID
  xyCoords <- inData[,c("LON","LAT")]
  sp.dat <- sp::SpatialPointsDataFrame(xyCoords, inData, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  sp.dat <- sp::spTransform(sp.dat, sp::CRS("+init=epsg:3857")) #epsg:3301
  return (as.data.frame(sp.dat))
}

getDbscan <- function( data, parameters ){
  allCoords <- dplyr::select(data, LON.1, LAT.1)
  eps <- calculateEps(allCoords,parameters$minPoints)
  return (dbscan::dbscan(allCoords, eps, parameters$minPoints))
}

calculateEps <- function( allCoords, minPoints ){
  distvec <- dbscan::kNNdist(allCoords, minPoints)
  distvec <- sort(distvec)
  distdf <- data.frame(c(1:length(distvec)),distvec)
  x <- distdf[,1]
  y <- distdf[,2]
  ede <- inflection::ede(x,y,0)
  return (y[ede[1]])
}

bindClusterToData <- function( dbscan, data ){
  cluster <- dbscan$cluster
  return (cbind(data, cluster, deparse.level = 1 ))
}

#apply PCA, return original data + flgagAll true/false
applyPCA <- function( data, parameters ){
#TODO make columns dynamic
  variables <- dplyr::select(data,HEIGHT,HEIGHT.WRT.DEM,SIGMA.HEIGHT,VEL,SIGMA.VEL,CUM.DISP,COHER)
  varMatrix <- as.matrix(variables)
  aR = rospca::robpca(varMatrix, k=2, ndir=5000)
  isCore = aR[["flag.all"]]
  return(cbind(data, isCore, deparse.level = 1))
}

getOCFreeData <- function( data ){
  return ( data[data$isCore == TRUE & data$cluster != 0,] )
}
#Get all pca which are not noise in dbscan
removeNoisePoints <- function( data ){
  return (data[data$isCore == FALSE & data$cluster != 0,])
}
getNoiseWithCoherLimit <- function( data, parameters ){
  return ( data[data$isCore == FALSE & data$cluster == 0 & data$COHER>parameters$minCoher,] )
}

#Application-driven outliers
#IN:
#data: pointsWithOCFlag
#eps: same as dbscan
#minCoher: default 0.7
#rejCrit:
#minJacc: index threshold default 0.6
processOCs <- function( pointsWithOCFlag, madOfVariablesPerClusterDF, parameters ){
  noiseFreeOutlierCandidates <- removeNoisePoints(pointsWithOCFlag)

  triang <- applyTriangulation( noiseFreeOutlierCandidates, parameters )
  delsgsOfGroupedOC <- triang$delsgs[triang$delsgs$dist<=parameters$eps,]
  groupedOCstoKeep <- processGroupedOutliers( noiseFreeOutlierCandidates, delsgsOfGroupedOC, madOfVariablesPerClusterDF, parameters )

  indsWithneighbors <- getIndsWithNeighbors(delsgsOfGroupedOC)
  indsOfIsolatedOutliers <- getIndsOfIsolatedOutliers( triang$ind.orig, indsWithneighbors ) #process isolated outliers

  isolatedOCstoKeep <- processIsolatedOutliers( noiseFreeOutlierCandidates, indsOfIsolatedOutliers, madOfVariablesPerClusterDF, parameters )

  return (rbind(groupedOCstoKeep,isolatedOCstoKeep))
  }

processIsolatedOutliers <- function( noiseFreeOutlierCandidates, indsOfIsolatedOutliers, madOfVariablesPerClusterDF, parameters ){
  idsToKeep <- vector()
  for (i in indsOfIsolatedOutliers){
    ocData <- noiseFreeOutlierCandidates[i,]
    if ( ocData$COHER>parameters$minCoher ){
      rejectionTable <- calculateMatrixPerOc(noiseFreeOutlierCandidates, i, madOfVariablesPerClusterDF)
      if(rowSums(rejectionTable) < 2){
        idsToKeep <- c(idsToKeep,i)
      }
    }
  }
  return( noiseFreeOutlierCandidates[idsToKeep,])
}

#####TODO lon1,lat1
#returns df with neighboring points and their distance
#TODO errors
applyTriangulation <- function( data, parameters ){
  xy <- data[,c("LON","LAT")]
  spdf <- sp::SpatialPointsDataFrame(coords = xy, data = data,
                                 proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) #
  # Voronoi
  triang <- deldir::deldir(data$LON.1, data$LAT.1 , digits=18) #aka triangulation
  triang$delsgs$dist <- raster::pointDistance(cbind(triang$delsgs$x1, triang$delsgs$y1),cbind(triang$delsgs$x2, triang$delsgs$y2), lonlat=FALSE)

  return(triang)
}

#Returns inds wchich have neighbors<=eps
getIndsWithNeighbors <- function(delsgsOfGroupedOC){
  return(sort.default(dplyr::union(delsgsOfGroupedOC$ind1,delsgsOfGroupedOC$ind2)))
}

getIndsOfGroupedOutliers <- function( indOrig, indsWithneighbors ){
  return (indOrig[(indOrig %in% indsWithneighbors)])
}
getIndsOfIsolatedOutliers <- function( indOrig, indsWithneighbors ){
  return ( indOrig[!(indOrig %in% indsWithneighbors)] )
}

processGroupedOutliers <- function( noiseFreeOutlierCandidates, delsgsOfGroupedOC, madOfVariablesPerClusterDF, parameters ){
  groupedComponents <- groupOCfromPCA( delsgsOfGroupedOC )
  return ( processEachGroup( noiseFreeOutlierCandidates, groupedComponents, madOfVariablesPerClusterDF, parameters ) )
}

#Makes graph of nonisolated OCs and groupes them. Returns groups with indices of original data.
#IN: allData <- points with oc and cluster flag. delsgsOfGroupedOC - delsgs only for grouped
groupOCfromPCA <- function( delsgsOfGroupedOC ){ #ocFromPcaWihtoutNoise
  groupedOC <- dplyr::select(delsgsOfGroupedOC, ind1, ind2, x1, y1, x2, y2, dist)
  g <- igraph::graph_from_data_frame(groupedOC, directed=FALSE )
  return ( igraph::groups( igraph::components( g ) ) )
}

processEachGroup <- function( noiseFreeOutlierCandidates, groupedComponents, madOfVariablesPerClusterDF, parameters ){
  pointsToKeep <- data.frame()
  for (i in 1:length(groupedComponents)){
    rejtable <-calculateMatrixPerOc(noiseFreeOutlierCandidates, groupedComponents[[i]], madOfVariablesPerClusterDF )
    pointsToKeep <- rbind(pointsToKeep,calculateJaccard(rejtable, noiseFreeOutlierCandidates, parameters))
  }
  return(pointsToKeep)
}


calculateMadOfVariablePerCluster <- function( pointsWithOCFlag, dbscan, parameters ){
  clustersWithMajOutliers <- getClustersWithMajorityOutliers( pointsWithOCFlag )
  groupCluster <- as.integer(names(which.max(table(dbscan$cluster))))
  cluster <- groupCluster
  clusters <-  as.integer(names(table(dbscan$cluster)))

  groupClusterData <- pointsWithOCFlag[pointsWithOCFlag$cluster == groupCluster & pointsWithOCFlag$isCore == TRUE,
                                       c("HEIGHT","HEIGHT.WRT.DEM","SIGMA.HEIGHT","VEL","SIGMA.VEL","CUM.DISP","COHER")]
  madPerColumnForGroupCluster <- as.data.frame(lapply(groupClusterData, madPerColumnFun, rj=parameters$rejCrit))
  result <- cbind(cluster, madPerColumnForGroupCluster )
  for (i in clusters[-1]){ #TODO selle saaks ilusamaks
    if ( i != groupCluster){
      cluster <- i
      if( cluster %in% clustersWithMajOutliers ){
        result <- rbind(result,cbind(cluster,madPerColumnForGroupCluster))
      } else {
        allNonOutliers <-  pointsWithOCFlag[pointsWithOCFlag$cluster == cluster & pointsWithOCFlag$isCore == TRUE,]
        allNonOutliers1 <-  allNonOutliers[,c("HEIGHT","HEIGHT.WRT.DEM","SIGMA.HEIGHT","VEL","SIGMA.VEL","CUM.DISP","COHER")]
        madPerColumn <-  as.data.frame(lapply(allNonOutliers1, madPerColumnFun, rj=parameters$rejCrit))
        result <- rbind(result,cbind(cluster,madPerColumn))
      }
    }
  }
  return(result)
}


calculateJaccard <- function(rejList, allData, parameters){
  ocGroupData <- allData[row.names(allData) %in% row.names(rejList),]
  jacc <- dist(rejList, method = "binary",
               diag = TRUE, upper = TRUE)
  jacc <- as.data.frame(as.matrix(jacc))
  jaccsim <- 1 - jacc
  jaccsim$id <- rownames(jaccsim)
  pairWiseJacc <- reshape2::melt(jaccsim, id.vars = "id")
  pairWiseJacc$id <- factor(pairWiseJacc$id, rownames(jaccsim))
  pairWiseJacc$variable <- factor(pairWiseJacc$variable, rev(rownames(jaccsim)))
  pairWiseJacc <- plyr::arrange(pairWiseJacc, variable, plyr::desc(id))
  different <- pairWiseJacc[pairWiseJacc$value < parameters$minJacc,]
  if (length(different[,1]) == 0){
    return(ocGroupData)
  }
  return(ocGroupData[row.names(ocGroupData) %in% unique(different$id) & ocGroupData$COHER> parameters$minCoher,])
}

getClustersWithMajorityOutliers <- function(data){
  library(dplyr)
  neg <- data[data$cluster!=0 & data$isCore==FALSE,]
  countNegByCluster <- neg %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(falseSum=n() )

  pos <- data[data$cluster!=0 & data$isCore==TRUE,]
  countPosByCluster <- pos %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(trueSum=n() )

  e <- dplyr::full_join(countPosByCluster,countNegByCluster, by="cluster")
  e <- e  %>% dplyr::mutate_all(funs(replace(., which(is.na(.)), 0)))

  return( e[e$falseSum > e$trueSum,]$cluster )
}

calculateMatrixPerOc <- function(noiseFreeOutlierCandidates, ocInds, madOfVariablesPerClusterDF ){ #viga, kui clusterid ei ole samad
  ocInds <- as.integer(ocInds)
  cluster <-  noiseFreeOutlierCandidates[ocInds[1],]$cluster

  madValues <- madOfVariablesPerClusterDF[madOfVariablesPerClusterDF$cluster == cluster, -1]

  rejectionTable <- noiseFreeOutlierCandidates[ocInds,c("HEIGHT","HEIGHT.WRT.DEM","SIGMA.HEIGHT","VEL","SIGMA.VEL","CUM.DISP","COHER")]
  for(i in 1:7){ #selle seitsme saab eraldi muutujaks length(väljad) vms
    madPerColumn <-
      rejectionTable[,i]<-ifelse(dplyr::between(rejectionTable[,i], madValues[1,i], madValues[2,i]),0,1)
  }
  return(rejectionTable)
}


madPerColumnFun <- function( x, rj ) {
  md <- mad(x,center=median(x),constant=1.4826)
  return(c(interval.lower(median(x),md,rj),interval.upper(median(x),md,rj)))
}
interval.upper <- function( median, mad, rj ) {
  median+rj*mad
}
interval.lower <- function( median, mad, rj ) {
  median-rj*mad
}
