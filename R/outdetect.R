extractOutliers <- function( data, minPoints, minCoher=0.7, k=2, rejCrit=3, minJacc = 0.6 ){
  params <- list(minPoints=minPoints, minCoher=minCoher, k=k, rejCrit=rejCrit, minJacc=minJacc)
  origColumnNames <- colnames(data)
  params$variables <- origColumnNames[4:length(colnames(data))]
  #validations

  allData <- prepareData(data)

  #Location based
  dbscan <- getDbscan( allData, params )
  params$eps <- dbscan$eps
  groupCluster <- as.integer(names(which.max(table(dbscan$cluster))))
  dataWithCluster <- bindClusterToData(dbscan, allData)

  #Data-driven
  pointsWithOCFlag <- applyPCA( dataWithCluster, params ) #points with oc flag; k no of components
  outlierAndNoiseFree <- getOCFreeData( pointsWithOCFlag ) #points to be included in final dataset and found by PCA
  noiseToBeIncluded <- getNoiseWithCoherLimit(pointsWithOCFlag, params)

  #Application-driven
  madOfVariablePerClusterDF <- calculateMadOfVariablePerCluster( pointsWithOCFlag, dbscan, params )
  additionalPoints <- processOCs( pointsWithOCFlag, madOfVariablePerClusterDF, params )
  pointsToKeep <- rbind(outlierAndNoiseFree, noiseToBeIncluded, additionalPoints)

  nonoutliers <- data[data$ID %in% row.names(pointsToKeep),]
  outliers <- data[!data$ID %in% row.names(pointsToKeep),]

  x <- list(nonoutliers = nonoutliers, outliers = outliers,params = params)

  return (x)
}


#TODO column + data validations
prepareData <- function( data ){
  inData <- data[-1]
  row.names(inData) <- data$ID
  xyCoords <- inData[,c("LON","LAT")]
  sp.dat <- sp::SpatialPointsDataFrame(xyCoords, inData, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  sp.dat <- sp::spTransform(sp.dat, sp::CRS("+init=epsg:3395"))
  return (as.data.frame(sp.dat))
}

getDbscan <- function( data, params ){
  allCoords <- data[,c("LON.1", "LAT.1")]
  eps <- calculateEps(allCoords,params$minPoints)
  return (dbscan::dbscan(allCoords, eps, params$minPoints))
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
applyPCA <- function( data, params ){
  variables <- data[,params$variables]

  varMatrix <- as.matrix(variables)
  aR = rrcov::PcaHubert(x=varMatrix, k=params$k, mcd=FALSE)
  isCore = aR@flag
  return(cbind(data, isCore, deparse.level = 1))
}

getOCFreeData <- function( data ){
  return ( data[data$isCore == TRUE, ] )
}
#Get all pca which are not noise in dbscan
removeNoisePoints <- function( data ){
  return (data[data$isCore == FALSE & data$cluster != 0,])
}
getNoiseWithCoherLimit <- function( data, params ){
  return ( data[data$isCore == FALSE & data$cluster == 0 & data$COHER>params$minCoher,] )
}

#Application-driven outliers
#IN:
#data: pointsWithOCFlag
#eps: same as dbscan
#minCoher: default 0.7
#rejCrit:
#minJacc: index threshold default 0.6
processOCs <- function( pointsWithOCFlag, madOfVariablesPerClusterDF, params ){
  noiseFreeOutlierCandidates <- removeNoisePoints(pointsWithOCFlag)
  triang <- applyTriangulation( noiseFreeOutlierCandidates, params )
  delsgsOfGroupedOC <- triang$delsgs[triang$delsgs$dist<=params$eps,]
  groupedOCstoKeep <- processGroupedOutliers( noiseFreeOutlierCandidates, delsgsOfGroupedOC, madOfVariablesPerClusterDF, params )

  indsWithneighbors <- getIndsWithNeighbors(delsgsOfGroupedOC)

  indsOfIsolatedOutliers <- getIndsOfIsolatedOutliers( triang$ind.orig, indsWithneighbors ) #process isolated outliers

  isolatedOCstoKeep <- processIsolatedOutliers( noiseFreeOutlierCandidates, indsOfIsolatedOutliers, madOfVariablesPerClusterDF, params )
  return (rbind(groupedOCstoKeep,isolatedOCstoKeep))
  }

processIsolatedOutliers <- function( noiseFreeOutlierCandidates, indsOfIsolatedOutliers, madOfVariablesPerClusterDF, params ){
  idsToKeep <- vector()
  for (i in indsOfIsolatedOutliers){
    ocData <- noiseFreeOutlierCandidates[i,]
    if ( ocData$COHER>params$minCoher ){
      rejectionTable <- calculateMatrixPerOc(noiseFreeOutlierCandidates, i, madOfVariablesPerClusterDF, params)
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
applyTriangulation <- function( data, params ){
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

processGroupedOutliers <- function( noiseFreeOutlierCandidates, delsgsOfGroupedOC, madOfVariablesPerClusterDF, params ){
  groupedComponents <- groupOCfromPCA( delsgsOfGroupedOC )
  return ( processEachGroup( noiseFreeOutlierCandidates, groupedComponents, madOfVariablesPerClusterDF, params ) )
}

#Makes graph of nonisolated OCs and groupes them. Returns groups with indices of original data.
#IN: allData <- points with oc and cluster flag. delsgsOfGroupedOC - delsgs only for grouped
groupOCfromPCA <- function( delsgsOfGroupedOC ){ #ocFromPcaWihtoutNoise
  groupedOC <- dplyr::select(delsgsOfGroupedOC, ind1, ind2, x1, y1, x2, y2, dist)
  g <- igraph::graph_from_data_frame(groupedOC, directed=FALSE )
  return ( igraph::groups( igraph::components( g ) ) )
}

processEachGroup <- function( noiseFreeOutlierCandidates, groupedComponents, madOfVariablesPerClusterDF, params ){
  pointsToKeep <- data.frame()
  for (i in 1:length(groupedComponents)){
    rejtable <-calculateMatrixPerOc(noiseFreeOutlierCandidates, groupedComponents[[i]], madOfVariablesPerClusterDF, params )
    pointsToKeep <- rbind(pointsToKeep,calculateJaccard(rejtable, noiseFreeOutlierCandidates, params))
  }
  return(pointsToKeep)
}


calculateMadOfVariablePerCluster <- function( pointsWithOCFlag, dbscan, params ){
  clustersWithMajOutliers <- getClustersWithMajorityOutliers( pointsWithOCFlag )
  groupCluster <- as.integer(names(which.max(table(dbscan$cluster))))
  cluster <- groupCluster
  clusters <-  as.integer(names(table(dbscan$cluster)))
  groupClusterData <- pointsWithOCFlag[pointsWithOCFlag$cluster == groupCluster & pointsWithOCFlag$isCore == TRUE,params$variables]
  madPerColumnForGroupCluster <- as.data.frame(lapply(groupClusterData, madPerColumnFun, rj=params$rejCrit))
  result <- cbind(cluster, madPerColumnForGroupCluster )
  for (i in clusters[-1]){ 
    if ( i != groupCluster){
      cluster <- i
      if( cluster %in% clustersWithMajOutliers ){
        result <- rbind(result,cbind(cluster,madPerColumnForGroupCluster))
      } else {
        allNonOutliers <-  pointsWithOCFlag[pointsWithOCFlag$cluster == cluster & pointsWithOCFlag$isCore == TRUE,]
        allNonOutliers1 <-  allNonOutliers[,params$variables]
        madPerColumn <-  as.data.frame(lapply(allNonOutliers1, madPerColumnFun, rj=params$rejCrit))
        result <- rbind(result,cbind(cluster,madPerColumn))
      }
    }
  }
  return(result)
}


calculateJaccard <- function(rejList, allData, params){
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
  different <- pairWiseJacc[pairWiseJacc$value < params$minJacc,]
  if (length(different[,1]) == 0){
    return(ocGroupData)
  }
  return(ocGroupData[row.names(ocGroupData) %in% unique(different$id) & ocGroupData$COHER> params$minCoher,])
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

calculateMatrixPerOc <- function(noiseFreeOutlierCandidates, ocInds, madOfVariablesPerClusterDF, params ){ #viga, kui clusterid ei ole samad
  ocInds <- as.integer(ocInds)
  cluster <-  noiseFreeOutlierCandidates[ocInds[1],]$cluster
  madValues <- madOfVariablesPerClusterDF[madOfVariablesPerClusterDF$cluster == cluster, -1]

  rejectionTable <- noiseFreeOutlierCandidates[ocInds,params$variables]
  for(i in 1:length(params$variables)){
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
