################ load necessary libraries ############################################################
library(igraph)
library(rNMF)
library(xgboost)

################ read input data #####################################################################
# read known miRNA-SM association dataset
knownMSA <- read.table(file = "./SM-miRNA/similar/SM-miRNA_Num_A_similar.csv", header = F,sep=",")

# # read miRNA functional similarity matrix
similaritiesOfMiRNA <- as.matrix(read.table(file = "./SM-miRNA/similar/miRNA_smilarity_maritx.csv", header = F,sep=","))
# # read SM similarity matrix
similaritiesOfSM <- as.matrix(read.table(file = "./SM-miRNA/similar/SM_similarity_matrix.csv", header = F,sep=","))



############### function to build training and testing data ##########################################

BuildTrainingAndTestingData <- function(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
                                        negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices, localLoocvTestingIndices) {
  
  ##############################
  ## Type 1 feature of miRNAs ##
  ##############################
  
  # number of observations in each row of MSA
  noOfObervationsOfMiRNA <- matrix(rep(0, m), nrow = m, ncol = 1)
    for(i in 1 : m) {
    noOfObervationsOfMiRNA[i] <- sum(MSA[i, ])
  }
  
  # average of all similarity scores for each miRNA
  aveOfSimilaritiesOfMiRNA <- matrix(rep(0, m), nrow = m, ncol = 1)
  for(i in 1 : m) {
    aveOfSimilaritiesOfMiRNA[i] <- mean(similaritiesOfMiRNA[i, ])
  }
  
  # histogram feature: cut [0, 1] into five bins and count the proportion of similarity scores that fall into each bin
  hist1MiRNA <- matrix(rep(0, m), nrow = m, ncol = 1) # [0, 0.2)
  hist2MiRNA <- matrix(rep(0, m), nrow = m, ncol = 1) # [0.2, 0.4)
  hist3MiRNA <- matrix(rep(0, m), nrow = m, ncol = 1) # [0.4, 0.6)
  hist4MiRNA <- matrix(rep(0, m), nrow = m, ncol = 1) # [0.6, 0.8)
  hist5MiRNA <- matrix(rep(0, m), nrow = m, ncol = 1) # [0.8, 1]
  for(i in 1: m) {
    hist1Count = 0
    hist2Count = 0
    hist3Count = 0
    hist4Count = 0
    hist5Count = 0
    for(j in 1 : m) {
      if(similaritiesOfMiRNA[i, j] < 0.2) {
        hist1Count = hist1Count + 1
      } else if(similaritiesOfMiRNA[i, j] < 0.4) {
        hist2Count = hist2Count + 1
      } else if(similaritiesOfMiRNA[i, j] < 0.6) {
        hist3Count = hist3Count + 1
      } else if(similaritiesOfMiRNA[i, j] < 0.8) {
        hist4Count = hist4Count + 1
      } else if(similaritiesOfMiRNA[i, j] <= 1) {
        hist5Count = hist5Count + 1
      }
    }
    hist1MiRNA[i] <- hist1Count / m
    hist2MiRNA[i] <- hist2Count / m
    hist3MiRNA[i] <- hist3Count / m
    hist4MiRNA[i] <- hist4Count / m
    hist5MiRNA[i] <- hist5Count / m
  }
  
  # concatenation
  feature1OfMiRNA <- cbind(noOfObervationsOfMiRNA, aveOfSimilaritiesOfMiRNA, hist1MiRNA, 
                           hist2MiRNA, hist3MiRNA, hist4MiRNA, hist5MiRNA)
  colnames(feature1OfMiRNA) <- c("noOfObervationsOfMiRNA", "aveOfSimilaritiesOfMiRNA", "hist1MiRNA", 
                                 "hist2MiRNA", "hist3MiRNA", "hist4MiRNA", 
                                 "hist5MiRNA")
  
  
  ################################
  ## Type 1 feature of SMs ##
  ################################
  
  # number of observations in each column of MSA
  noOfObervationsOfSM <- matrix(rep(0, s), nrow = s, ncol = 1)
  for(i in 1 : s) {
    noOfObervationsOfSM[i] <- sum(MSA[, i])
  }
  
  # average of all similarity scores for each SM
  aveOfSimilaritiesOfSM <- matrix(rep(0, s), nrow = s, ncol = 1)
  for(i in 1 : s) {
    aveOfSimilaritiesOfSM[i] <- mean(similaritiesOfSM[, i])
  }
  
  # histogram feature: cut [0, 1] into five bins and count the proportion of similarity scores that fall into each bin
  hist1SM <- matrix(rep(0, s), nrow = s, ncol = 1) # [0, 0.2)
  hist2SM <- matrix(rep(0, s), nrow = s, ncol = 1) # [0.2, 0.4)
  hist3SM <- matrix(rep(0, s), nrow = s, ncol = 1) # [0.4, 0.6)
  hist4SM <- matrix(rep(0, s), nrow = s, ncol = 1) # [0.6, 0.8)
  hist5SM <- matrix(rep(0, s), nrow = s, ncol = 1) # [0.8, 1]
  for(i in 1: s) {
    hist1Count = 0
    hist2Count = 0
    hist3Count = 0
    hist4Count = 0
    hist5Count = 0
    for(j in 1 : s) {
      if(similaritiesOfSM[i, j] < 0.2) {
        hist1Count = hist1Count + 1
      } else if(similaritiesOfSM[i, j] < 0.4) {
        hist2Count = hist2Count + 1
      } else if(similaritiesOfSM[i, j] < 0.6) {
        hist3Count = hist3Count + 1
      } else if(similaritiesOfSM[i, j] < 0.8) {
        hist4Count = hist4Count + 1
      } else if(similaritiesOfSM[i, j] <= 1) {
        hist5Count = hist5Count + 1
      }
    }
    hist1SM[i] <- hist1Count / s
    hist2SM[i] <- hist2Count / s
    hist3SM[i] <- hist3Count / s
    hist4SM[i] <- hist4Count / s
    hist5SM[i] <- hist5Count / s
  }
  
  # concatenation 
  feature1OfSM <- cbind(noOfObervationsOfSM, aveOfSimilaritiesOfSM, hist1SM, 
                             hist2SM, hist3SM, hist4SM, hist5SM)
  colnames(feature1OfSM) <- c("noOfObervationsOfSM", "aveOfSimilaritiesOfSM", "hist1SM", 
                                   "hist2SM", "hist3SM", "hist4SM", 
                                   "hist5SM")
  
  
  ##############################
  ## Type 2 feature of miRNAs ##
  ##############################
  
  # number of neighbors of miRNAs and similarity values for 10 nearest neighbors
  numberOfNeighborsMiRNA <- matrix(rep(0, m), nrow = m, ncol = 1)
  similarities10KnnMiRNA <- matrix(rep(0, 10 * m), nrow = m, ncol = 10)
  averageOfFeature1MiRNA <- matrix(rep(0, 7 * m), nrow = m, ncol = 7)
  weightedAverageOfFeature1MiRNA <- matrix(rep(0, 7 * m), nrow = m, ncol = 7)
  similarityGraphMiRNA <- matrix(rep(0, m * m), nrow = m, ncol = m)
  meanSimilarityMiRNA <- mean(similaritiesOfMiRNA)
  for(i in 1 : m) {
    neighborCount = 0 - 1 # similarity between an miRNA and itself is not counted 
    for(j in 1 : m) {
      if(similaritiesOfMiRNA[i, j] >= meanSimilarityMiRNA) {
        neighborCount = neighborCount + 1
        similarityGraphMiRNA[i, j] = 1
      }
    }
    numberOfNeighborsMiRNA[i] <- neighborCount
    similarities10KnnMiRNA[i, ] <- sort(similaritiesOfMiRNA[i, ], decreasing = T, index.return = T)$x[2:11]
    indices <- sort(similaritiesOfMiRNA[i, ], decreasing = T, index.return = T)$ix[2:11]
    if(neighborCount == 0) {
      averageOfFeature1MiRNA[i, ] <- rep(0, 7)
      weightedAverageOfFeature1MiRNA[i, ] <- rep(0, 7)
      next
    } else if(neighborCount == 1) {
      averageOfFeature1MiRNA[i, ] <- feature1OfMiRNA[indices[1], ] / 10
      weightedAverageOfFeature1MiRNA[i, ] <- feature1OfMiRNA[indices[1], ] * similarities10KnnMiRNA[i, ][1] / 10
      next
    } else if (neighborCount <= length(indices)) {
      indices <- indices[1 : neighborCount]
    }
    averageOfFeature1MiRNA[i, ] <- apply(feature1OfMiRNA[indices, ], MARGIN = 2, FUN = function(x) sum(x) / 10) # divide by 10 to make the mean calculation fair for those miRNAs with less than 10 neighbors
    weightedAverageOfFeature1MiRNA[i, ] <- apply(feature1OfMiRNA[indices, ], MARGIN = 2, 
                                                 FUN = function(x) sum(x * similarities10KnnMiRNA[i, ][1 : length(indices)]) / 10)
  }
  
  # build miRNA similarity graph
  similarityIgraphMiRNA <- graph_from_adjacency_matrix(adjmatrix = similarityGraphMiRNA, mode = "undirected", weighted = NULL,
                                                       diag = T)
  betweennessCentralityMiRNA <- betweenness(similarityIgraphMiRNA, directed = F, normalized = T)
  closenessCentralityMiRNA <- closeness(similarityIgraphMiRNA, mode = "all")
  eigenVectorCentralityMiRNA <- eigen_centrality(similarityIgraphMiRNA, directed = F)$vector
  pageRankMiRNA <- page.rank(similarityIgraphMiRNA, directed = F)$vector

  # concatenation
  feature2OfMiRNA <- cbind(numberOfNeighborsMiRNA, similarities10KnnMiRNA, averageOfFeature1MiRNA, weightedAverageOfFeature1MiRNA,
                           betweennessCentralityMiRNA, closenessCentralityMiRNA, eigenVectorCentralityMiRNA, pageRankMiRNA)
  colnames(feature2OfMiRNA) <- c("numberOfNeighborsMiRNA", "knn1SimilarityMiRNA", "knn2SimilarityMiRNA", "knn3SimilarityMiRNA", 
                                 "knn4SimilarityMiRNA", "knn5SimilarityMiRNA", "knn6SimilarityMiRNA", "knn7SimilarityMiRNA",
                                 "knn8SimilarityMiRNA", "knn9SimilarityMiRNA", "knn10SimilarityMiRNA", "aveNoObsMiRNA", 
                                 "aveOfAveSimilarityMiRNA", "aveHist1MiRNA", "aveHist2MiRNA", "aveHist3MiRNA", "aveHist4MiRNA", 
                                 "aveHist5MiRNA", "weightedAveNoObsMiRNA", "weightedAveOfAveSimilarityMiRNA", "weightedAveHist1MiRNA",
                                 "weightedAveHist2MiRNA", "weightedAveHist3MiRNA", "weightedAveHist4MiRNA", "weightedAveHist5MiRNA", 
                                 "betweennessCentralityMiRNA", "closenessCentralityMiRNA", "eigenVectorCentralityMiRNA", "pageRankMiRNA")
  
  
  ################################
  ## Type 2 feature of SMs ##
  ################################
  
  # number of neighbors of SMs and similarity values for 10 nearest neighbors
  numberOfNeighborsSM <- matrix(rep(0, s), nrow = s, ncol = 1)
  similarities10KnnSM <- matrix(rep(0, 10 * s), nrow = s, ncol = 10)
  averageOfFeature1SM <- matrix(rep(0, 7 * s), nrow = s, ncol = 7)
  weightedAverageOfFeature1SM <- matrix(rep(0, 7 * s), nrow = s, ncol = 7)
  similarityGraphSM <- matrix(rep(0, s * s), nrow = s, ncol = s)
  meanSimilaritySM <- mean(similaritiesOfSM)
  for(i in 1 : s) {
    neighborCount = 0 - 1 # similarity between a SM and itself is not counted 
    for(j in 1 : s) {
      if(similaritiesOfSM[i, j] >= meanSimilaritySM) {
        neighborCount = neighborCount + 1
        similarityGraphSM[i, j] = 1
      }
    }
    numberOfNeighborsSM[i] <- neighborCount
    similarities10KnnSM[i, ] <- sort(similaritiesOfSM[i, ], decreasing = T, index.return = T)$x[2:11]
    indices <- sort(similaritiesOfSM[i, ], decreasing = T, index.return = T)$ix[2:11]
    if(neighborCount == 0) {
      averageOfFeature1SM[i, ] <- rep(0, 7)
      weightedAverageOfFeature1SM[i, ] <- rep(0, 7)
      next
    } else if(neighborCount == 1) {
      averageOfFeature1SM[i, ] <- feature1OfSM[indices[1], ] / 10
      weightedAverageOfFeature1SM[i, ] <- feature1OfSM[indices[1], ] * similarities10KnnSM[i, ][1] / 10
      next
    } else if (neighborCount <= length(indices)) {
      indices <- indices[1 : neighborCount]
    }
    averageOfFeature1SM[i, ] <- apply(feature1OfSM[indices, ], MARGIN = 2, FUN = function(x) sum(x) / 10) # divide by 10 to make the mean calculation fair for those SMs with less than 10 neighbors
    weightedAverageOfFeature1SM[i, ] <- apply(feature1OfSM[indices, ], MARGIN = 2, 
                                                   FUN = function(x) sum(x * similarities10KnnSM[i, ][1 : length(indices)]) / 10)
  }
  
  # build SM similarity graph
  library(igraph)
  similarityIgraphSM <- graph_from_adjacency_matrix(adjmatrix = similarityGraphSM, mode = "undirected", weighted = NULL,
                                                         diag = T)
  betweennessCentralitySM <- betweenness(similarityIgraphSM, directed = F, normalized = T)
  closenessCentralitySM <- closeness(similarityIgraphSM, mode = "all")
  eigenVectorCentralitySM <- eigen_centrality(similarityIgraphSM, directed = F)$vector
  pageRankSM <- page.rank(similarityIgraphSM, directed = F)$vector

  # concatenation
  feature2OfSM <- cbind(numberOfNeighborsSM, similarities10KnnSM, averageOfFeature1SM, weightedAverageOfFeature1SM,
                             betweennessCentralitySM, closenessCentralitySM, eigenVectorCentralitySM, pageRankSM)
  colnames(feature2OfSM) <- c("numberOfNeighborsSM", "knn1SimilaritySM", "knn2SimilaritySM", "knn3SimilaritySM", 
                                   "knn4SimilaritySM", "knn5SimilaritySM", "knn6SimilaritySM", "knn7SimilaritySM",
                                   "knn8SimilaritySM", "knn9SimilaritySM", "knn10SimilaritySM", "aveNoObsSM", 
                                   "aveOfAveSimilaritySM", "aveHist1SM", "aveHist2SM", "aveHist3SM", "aveHist4SM", 
                                   "aveHist5SM", "weightedAveNoObsSM", "weightedAveOfAveSimilaritySM", "weightedAveHist1SM",
                                   "weightedAveHist2SM", "weightedAveHist3SM", "weightedAveHist4SM", "weightedAveHist5SM", 
                                   "betweennessCentralitySM", "closenessCentralitySM", "eigenVectorCentralitySM", "pageRankSM")
  
  
  ###########################################
  ## Type 3 feature of miRNA-SM pairs ##
  ###########################################
  
  # matrix factorization
  set.seed(666)
  mfMSA <- rnmf(MSA, quiet = T, showprogress = F)
  latentVectorsMiRNA <- mfMSA$W
  latentVectorsSM <- mfMSA$H 
  # number of associations between an miRNA and a SM's neighbors
  numberOfSMNeighborAssociations <- c(rep(0, m * s))
  for(i in 1 : m) {
    for(j in 1 : s) {
      numberOfAssociations = ifelse(MSA[i, j] == 1, -1, 0)
      SMNeighbors = which(t(similaritiesOfSM[j, ]) >= meanSimilaritySM, arr.ind = F)
      for(k in 1 : length(SMNeighbors)) {
        if(MSA[i, SMNeighbors[k]] == 1) {
          numberOfAssociations = numberOfAssociations + 1
        }
      }
      numberOfSMNeighborAssociations[(i-1)*s+j] <- numberOfAssociations
    }
  }
  
  # number of associations between a SM and an miRNA's neighbors
  numberOfMiRNANeighborAssociations <- c(rep(0, m * s))
  for(i in 1 : s) {
    for(j in 1 : m) {
      numberOfAssociations = ifelse(MSA[j, i] == 1, -1, 0)
      miRNANeighbors = which(t(similaritiesOfMiRNA[j, ]) >= meanSimilarityMiRNA, arr.ind = F)
      for(k in 1 : length(miRNANeighbors)) {
        if(MSA[miRNANeighbors[k]] == 1) {
          numberOfAssociations = numberOfAssociations + 1
        }
      }
      numberOfMiRNANeighborAssociations[(i-1)*m+j] <- numberOfAssociations
    }
  }
  
  # build MSA graph
  MSAGraph <- graph_from_incidence_matrix(incidence = MSA, directed = F, mode = "total")
  betweennessCentralityMSA <- betweenness(MSAGraph, directed = F, normalized = T)
  betweennessCentralityMiRNAInMSA <- betweennessCentralityMSA[1:541]
  betweennessCentralitySMInMSA <- betweennessCentralityMSA[542:1372]
  closenessCentralityMSA <- closeness(MSAGraph, mode = "all")
  closenessCentralityMiRNAInMSA <- closenessCentralityMSA[1:541]
  closenessCentralitySMInMSA <- closenessCentralityMSA[542:1372]
  eigenVectorCentralityMSA <- eigen_centrality(MSAGraph, directed = F)$vector
  eigenVectorCentralityMiRNAInMSA <- eigenVectorCentralityMSA[1:541]
  eigenVectorCentralitySMInMSA <- eigenVectorCentralityMSA[542:1372]
  pageRankMSA <- page.rank(MSAGraph, directed = F)$vector
  pageRankMiRNAInMSA <- pageRankMSA[1:541]
  pageRankSMInMSA <- pageRankMSA[542:1372]
 
  #########################################
  ## function to combine feature vectors ##
  #########################################
  
  BuildFeatures <- function(positiveAndNegativeIndices) {
    positiveAndNegativeMiRNAIndices <- ifelse(positiveAndNegativeIndices %% 831 == 0, positiveAndNegativeIndices / 831, 
                                              as.integer(positiveAndNegativeIndices / 831) + 1) 
    positiveAndNegativeSMIndices <- ifelse(positiveAndNegativeIndices %% 831 == 0, 831, positiveAndNegativeIndices %% 831)
    loocvFeature1MiRNA <- feature1OfMiRNA[positiveAndNegativeMiRNAIndices, ]
    loocvFeature2MiRNA <- feature2OfMiRNA[positiveAndNegativeMiRNAIndices, ]
    loocvFeature1SM <- feature1OfSM[positiveAndNegativeSMIndices, ]
    loocvFeature2SM <- feature2OfSM[positiveAndNegativeSMIndices, ]
    loocvFeature3 <- cbind(latentVectorsMiRNA[positiveAndNegativeMiRNAIndices, ], t(latentVectorsSM[, positiveAndNegativeSMIndices]), 
                           numberOfSMNeighborAssociations[positiveAndNegativeMiRNAIndices], numberOfMiRNANeighborAssociations[positiveAndNegativeSMIndices], 
                           betweennessCentralityMiRNAInMSA[positiveAndNegativeMiRNAIndices], closenessCentralityMiRNAInMSA[positiveAndNegativeMiRNAIndices], 
                           eigenVectorCentralityMiRNAInMSA[positiveAndNegativeMiRNAIndices], pageRankMiRNAInMSA[positiveAndNegativeMiRNAIndices], 
                           betweennessCentralitySMInMSA[positiveAndNegativeSMIndices], closenessCentralitySMInMSA[positiveAndNegativeSMIndices],
                           eigenVectorCentralitySMInMSA[positiveAndNegativeSMIndices], pageRankSMInMSA[positiveAndNegativeSMIndices])
    colnames(loocvFeature3) <- c("latentVectors1MiRNA", "latentVectors2MiRNA", "latentVectors3MiRNA", "latentVectors4MiRNA", "latentVectors5MiRNA", "latentVectors1SM",
                                 "latentVectors2SM", "latentVectors3SM", "latentVectors4SM", "latentVectors5SM", "numberOfSMNeighborAssociations", 
                                 "numberOfMiRNANeighborAssociations", "betweennessCentralityMiRNAInMSA", "closenessCentralityMiRNAInMSA", "eigenVectorCentralityMiRNAInMSA", 
                                 "pageRankMiRNAInMSA", "betweennessCentralitySMInMSA", "closenessCentralitySMInMSA", "eigenVectorCentralitySMInMSA",
                                 "pageRankSMInMSA")
    loocvFeatureVectors <- cbind(loocvFeature1MiRNA, loocvFeature1SM, loocvFeature2MiRNA, loocvFeature2SM, loocvFeature3)
    return(loocvFeatureVectors)
  }
  
  # build training labels
  trainingLabels <- matrix(c(rep(1, length(knownMSAIndices)), rep(0, length(negativeSampleIndices))), nrow = length(knownMSAIndices) + 
                             length(negativeSampleIndices), ncol = 1)
  colnames(trainingLabels) <- "labels"
  
  # build loocv training data
  loocvTrainingFeatureVectors <- cbind(trainingLabels, BuildFeatures(positiveAndNegativeIndices))
  
  # build global loocv testing data
  globalLoocvTestingFeatureVectors <- BuildFeatures(globalLoocvTestingIndices)
  
  # build local loocv testing data
  localLoocvTestingFeatureVectors <- BuildFeatures(localLoocvTestingIndices)
  
  return(list("loocvTrainingFeatureVectors" = loocvTrainingFeatureVectors, 
              "globalLoocvTestingFeatureVectors" = globalLoocvTestingFeatureVectors, 
              "localLoocvTestingFeatureVectors" = localLoocvTestingFeatureVectors))
}

############### loops for global and local loocv #####################################################
# record the miRNA count, SM count and known MSA count
m <- 541 # 495
s <- 831 # 383
noOfKnownMSA <- nrow(knownMSA) # 664

# selection of xgboost parameters
parameters <- list(eta = 1, maxDepth = 6, lambda = 1, gamma = 0)

# placeholder for loocv rankings
rankings <- matrix(nrow = 0, ncol = 2)
colnames(rankings) <- c("globalRankings", "localRankings")

# loocv loops
for(negated in 1 : nrow(knownMSA)) {
  
  # find negated miRNA, SM and their association's index
  negatedMiRNA <- knownMSA$V2[negated]
  negatedSM <- knownMSA$V1[negated]
  negatedIndex <- (negatedMiRNA - 1) * s + negatedSM
  
  # build MSA matrix
  loocvKnownMSA <- knownMSA[-negated, ]
  originalMSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  for(i in 1 : m) {
    negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
    for(j in 1 : s) {
      if (j %in% negatedAssociations$V1) {
        originalMSA[i, j] <-  1
      }
    }
  }

  # randomly select 663 negative samples
  knownMSAIndices <- which(t(originalMSA) == 1, arr.ind = F)
  allIndices <- 1 : (m * s)
  set.seed(666)
  negativeSampleIndices <- sample(allIndices[-knownMSAIndices], size = 663, replace = F)
  
  # find indices for training data
  positiveAndNegativeIndices <- c(knownMSAIndices, negativeSampleIndices)
  
  # find indices for global and local testing data
  globalLoocvTestingIndices <- (1 : (m * s))[-knownMSAIndices]
  negatedIndexInGlobalTesting <- which(globalLoocvTestingIndices == negatedIndex)
  negatedIndexInLocalTesting <-  which(which(originalMSA[negatedMiRNA,] == 0) == negatedSM)
  localLoocvTestingIndices <- (which(originalMSA[negatedMiRNA,] == 0) - 1) * m + negatedMiRNA
  #negatedIndexInLocalTesting <- which(which(MSA[,negatedSM] == 0) == negatedMiRNA)
  #localLoocvTestingIndices <- (which(MSA[,negatedSM] == 0) - 1) * s + negatedSM   

  #MSA preprocessing
  MSA<-originalMSA
  # for(i in 1:m){
  #   for(j in 1:s){
  #     if(originalMSA[i,j]==0){
  #       ##find similar MiRNA by LINEs
  #       SimilarityDeacreaseOfMiRNA=sort.int(similaritiesOfMiRNA[i,],decreasing=T,index.return=T)
  #       count1=0                    #count the top similar MiRNAs
  #       count2=0
  #       similarMiRNAIndex<-rep(0,3) #preset index of top3 similar miRNA
  #       similarMiRNASValue<-rep(0,3)#preset similarity value of top3 similar miRNA
  #       for(k in 1:m){                             #find similar MiRNA by LINEs
  #           # flag=ifelse(SimilarityDeacreaseOfMiRNA$x[k]>0.5,1,-1)
  #           # while(flag==1&&count1<3){               #only find top3 similar MiRNAs with similarity value above 0.5
  #             if(originalMSA[k,j]==1&&count1<3){      #find top3 similar SM in all other MiRNAs
  #              count1=count1+1
  #              similarMiRNAIndex[count1]=SimilarityDeacreaseOfMiRNA$ix[k]
  #              similarMiRNASValue[count1]=SimilarityDeacreaseOfMiRNA$x[k]
  #             }
  #             else next
  #         }
  #         predictSVsepMiRNA<-rep(0,3)
  #         for(l in 1:3){
  #           if(is.nan(similarMiRNASValue[l])){similarMiRNASValue[l]=0}
  #           predictSVsepMiRNA[l]=(similarMiRNASValue[l]^2)/sum(similarMiRNASValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
  #           if(is.nan(predictSVsepMiRNA[l])){predictSVsepMiRNA[l]=0}
  #           else next
  #         }
  #         predictSVMiRNA=sum(predictSVsepMiRNA)
  #         ##find similar SM by ROWs##
  #         SimilarityDeacreaseOfSM=sort.int(similaritiesOfSM[,j],decreasing=T,index.return=T)
  #         similarSMIndex<-rep(0,3) #preset index of top3 similar SM
  #         similarSMSValue<-rep(0,3)#preset similarity value of top3 similar SM
  #         for(k in 1:s){                           #find similar SM by ROWs
  #          # flag=ifelse(SimilarityDeacreaseOfSM$x[k]>0.5,1,-1)
  #          # while(flag==1&&count2<3){               #only find top3 similar SM with similarity value above 0.5                                
  #             if(originalMSA[i,k]==1&&count2<3){      #find top3 similar SM in all other SMs
  #              count2=count2+1
  #              similarSMIndex[count2]=SimilarityDeacreaseOfSM$ix[k]
  #              similarSMSValue[count2]=SimilarityDeacreaseOfSM$x[k]
  #             }
  #             else next
  #           }
  #   #}    
  #         predictSVsepSM<-rep(0,3)
  #         for(l in 1:3){
  #           if(is.nan(similarSMSValue[l])){similarSMSValue[l]=0}
  #           predictSVsepSM[l]=(similarSMSValue[l]^2)/sum(similarSMSValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
  #           if(is.nan(predictSVsepSM[l])){predictSVsepSM[l]=0}
  #           else next
  #         }
  #         predictSVSM=sum(predictSVsepSM)
  #         MSA[i,j]=(predictSVMiRNA+predictSVSM)/2
  #       }
        
  #       else next
  #       }
  #   }

##########################
  # MSA<-sMSA
  # for(i in 1:m){
  #   for(j in 1:s){
  #     if(sMSA[i,j]==0){
  #       ##find similar MiRNA by LINEs
  #       SimilarityDeacreaseOfMiRNA=sort.int(similaritiesOfMiRNA[i,],decreasing=T,index.return=T)
  #       count1=0                    #count the top similar MiRNAs
  #       count2=0
  #       similarMiRNAIndex<-rep(0,3) #preset index of top3 similar miRNA
  #       similarMiRNASValue<-rep(0,3)#preset similarity value of top3 similar miRNA
  #       for(k in 1:m){                             #find similar MiRNA by LINEs
  #           # flag=ifelse(SimilarityDeacreaseOfMiRNA$x[k]>0.5,1,-1)
  #           # while(flag==1&&count1<3){               #only find top3 similar MiRNAs with similarity value above 0.5
  #             if(sMSA[k,j]==1&&count1<3){      #find top3 similar SM in all other MiRNAs
  #              count1=count1+1
  #              similarMiRNAIndex[count1]=SimilarityDeacreaseOfMiRNA$ix[k]
  #              similarMiRNASValue[count1]=SimilarityDeacreaseOfMiRNA$x[k]
  #             }
  #             else next
  #         }
  #         predictSVsepMiRNA<-rep(0,3)
  #         for(l in 1:3){
  #           if(is.nan(similarMiRNASValue[l])){similarMiRNASValue[l]=0}
  #           predictSVsepMiRNA[l]=(similarMiRNASValue[l]^2)/sum(similarMiRNASValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
  #           if(is.nan(predictSVsepMiRNA[l])){predictSVsepMiRNA[l]=0}
  #           else next
  #         }
  #         predictSVMiRNA=sum(predictSVsepMiRNA)
  #         ##find similar SM by ROWs##
  #         SimilarityDeacreaseOfSM=sort.int(similaritiesOfSM[,j],decreasing=T,index.return=T)
  #         similarSMIndex<-rep(0,3) #preset index of top3 similar SM
  #         similarSMSValue<-rep(0,3)#preset similarity value of top3 similar SM
  #         for(k in 1:s){                           #find similar SM by ROWs
  #          # flag=ifelse(SimilarityDeacreaseOfSM$x[k]>0.5,1,-1)
  #          # while(flag==1&&count2<3){               #only find top3 similar SM with similarity value above 0.5                                
  #             if(sMSA[i,k]==1&&count2<3){      #find top3 similar SM in all other SMs
  #              count1=count1+1
  #              similarSMIndex[count2]=SimilarityDeacreaseOfSM$ix[k]
  #              similarSMSValue[count2]=SimilarityDeacreaseOfSM$x[k]
  #             }
  #             else next
  #           }
  #   #}    
  #         predictSVsepSM<-rep(0,3)
  #         for(l in 1:3){
  #           if(is.nan(similarSMSValue[l])){similarSMSValue[l]=0}
  #           predictSVsepSM[l]=(similarSMSValue[l]^2)/sum(similarSMSValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
  #           if(is.nan(predictSVsepSM[l])){predictSVsepSM[l]=0}
  #           else next
  #         }
  #         predictSVSM=sum(predictSVsepSM)
  #         MSA[i,j]=(predictSVMiRNA+predictSVSM)/2
  #       }
        
  #       else next
  #       }
  #   }
##########################

  # build training and testing data
  trainingAndTestingData <- BuildTrainingAndTestingData(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
                                                        negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices,
                                                        localLoocvTestingIndices)
  # fit xgboost
  xgboostLoocv <- xgboost(data = trainingAndTestingData$loocvTrainingFeatureVectors[,-1], booster = "gbtree", 
                          label = trainingAndTestingData$loocvTrainingFeatureVectors[,1], params = parameters, nthread = 2, nrounds = 1, 
                          objective = "binary:logitraw",verbose=2)
  
  # prediction
  predictedWeightsGlobal <- predict(xgboostLoocv, trainingAndTestingData$globalLoocvTestingFeatureVectors)
  predictedWeightsLocal <- predict(xgboostLoocv, trainingAndTestingData$localLoocvTestingFeatureVectors)
  
  # build rankings
  ########
  #average_ranking
  ########

  # globalscoreDicrease=sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)
  # localscoreDicrease=sort.int(predictedWeightsLocal, decreasing = T, index.return = T)

  # globalRankingOfNegated <- which(globalscoreDicrease$ix == negatedIndexInGlobalTesting)
  # localRankingOfNegated <- which(localscoreDicrease$ix == negatedIndexInLocalTesting)
  
  # globalfinalrank <- mean(which(globalscoreDicrease$x==globalscoreDicrease$x[globalRankingOfNegated]))
  # localfinalrank <- mean(which(localscoreDicrease$x==localscoreDicrease$x[localRankingOfNegated]))
 
  # rankings <- rbind(rankings, c(globalfinalrank, localfinalrank))

   ############
  #accurate ranking
  ############
  globalRankingOfNegated <- which(sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)$ix == negatedIndexInGlobalTesting)
  localRankingOfNegated <- which(sort.int(predictedWeightsLocal, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting)
  rankings <- rbind(rankings, c(globalRankingOfNegated, localRankingOfNegated))

}

# write rankings to disk
write.table(rankings, file = "./global_and_local_loocv_rankings.cvs", row.names = F)
