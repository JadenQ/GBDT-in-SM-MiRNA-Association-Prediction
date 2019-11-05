################ load necessary libraries ############################################################
library(igraph)
library(rNMF)
library(xgboost)
library(fpc)
library(cluster)
library(parallel)
library(pbapply)
################ read input data #####################################################################

#similar_group
####read known miRNA-SM association dataset
# knownMSA <- read.table(file = "./SM-miRNA/similar/SM-miRNA_Num_A_similar.csv", header = F,sep=",")

# # # read miRNA functional similarity matrix
# similaritiesOfMiRNA <- as.matrix(read.table(file = "./SM-miRNA/similar/miRNA_smilarity_maritx.csv", header = F,sep=","))
# # read SM similarity matrix
# similaritiesOfSM <- as.matrix(read.table(file = "./SM-miRNA/similar/SM_similarity_matrix.csv", header = F,sep=","))

# related_group
# read known miRNA-SM association dataset
knownMSA <- read.table(file = "./SM-miRNA/related/SM-miRNA_Num_A.csv", header = F,sep=",")

# read miRNA functional similarity matrix
similaritiesOfMiRNA <- as.matrix(read.table(file = "./SM-miRNA/related/miRNA_similarity_matrix2.csv", header = F,sep=","))

# read SM similarity matrix
similaritiesOfSM <- as.matrix(read.table(file = "./SM-miRNA/related/SM_similarity_matrix2.csv", header = F,sep=","))



MshapeOfMiRNA=dim(similaritiesOfMiRNA)

ALL_MiRNA=MshapeOfMiRNA[c(1)]

MshapeOfSM=dim(similaritiesOfSM)

ALL_SM=MshapeOfSM[c(1)]

# similaritiesOfMiRNA<-similaritiesOfMiRNA/0.7049
# similaritiesOfMiRNA[ALL_MiRNA,ALL_MiRNA]=similaritiesOfMiRNA[ALL_MiRNA,ALL_MiRNA]*0.7049
# similaritiesOfSM <- similaritiesOfSM/0.8906

############### function to build training and testing data ##########################################

BuildTrainingAndTestingData <- function(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
                                        negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices, localLoocvTestingIndices,localLoocvTestingIndices2) {
  
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
  closenessCentralitySM <- closeness(similarityIgraphSM, mode ="all")
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
  mfMSA <- rnmf(MSA,k=5,quiet = T, showprogress = F)
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
  betweennessCentralityMiRNAInMSA <- betweennessCentralityMSA[1:ALL_MiRNA]
  betweennessCentralitySMInMSA <- betweennessCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  closenessCentralityMSA <- closeness(MSAGraph, mode = "all")
  closenessCentralityMiRNAInMSA <- closenessCentralityMSA[1:ALL_MiRNA]
  closenessCentralitySMInMSA <- closenessCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  eigenVectorCentralityMSA <- eigen_centrality(MSAGraph, directed = F)$vector
  eigenVectorCentralityMiRNAInMSA <- eigenVectorCentralityMSA[1:ALL_MiRNA]
  eigenVectorCentralitySMInMSA <- eigenVectorCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  pageRankMSA <- page.rank(MSAGraph, directed = F)$vector
  pageRankMiRNAInMSA <- pageRankMSA[1:ALL_MiRNA]
  pageRankSMInMSA <- pageRankMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  #########################################
  ## function to combine feature vectors ##
  #########################################
  
  BuildFeatures <- function(positiveAndNegativeIndices) {
    positiveAndNegativeMiRNAIndices <- ifelse(positiveAndNegativeIndices %% ALL_SM == 0, positiveAndNegativeIndices / ALL_SM, 
                                              as.integer(positiveAndNegativeIndices / ALL_SM) + 1) 
    positiveAndNegativeSMIndices <- ifelse(positiveAndNegativeIndices %% ALL_SM == 0, ALL_SM, positiveAndNegativeIndices %% ALL_SM)
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
    # colnames(loocvFeature3) <- c("latentVectors1MiRNA", "latentVectors2MiRNA", "latentVectors3MiRNA", "latentVectors4MiRNA", "latentVectors5MiRNA","latentVectors6MiRNA", "latentVectors7MiRNA", "latentVectors8MiRNA", "latentVectors9MiRNA", "latentVectors10MiRNA",
    #                              "latentVectors11MiRNA", "latentVectors12MiRNA", "latentVectors13MiRNA", "latentVectors14MiRNA", "latentVectors15MiRNA","latentVectors16MiRNA", "latentVectors17MiRNA", "latentVectors18MiRNA", "latentVectors19MiRNA", "latentVectors20MiRNA",
    #                              "latentVectors21MiRNA", "latentVectors22MiRNA", "latentVectors23MiRNA", "latentVectors24MiRNA", "latentVectors25MiRNA","latentVectors26MiRNA", "latentVectors27MiRNA", "latentVectors28MiRNA", "latentVectors29MiRNA", "latentVectors30MiRNA",
    #                              "latentVectors1SM","latentVectors2SM", "latentVectors3SM", "latentVectors4SM", "latentVectors5SM","latentVectors6SM","latentVectors7SM", "latentVectors8SM", "latentVectors9SM", "latentVectors10SM",
    #                              "latentVectors11SM","latentVectors12SM", "latentVectors13SM", "latentVectors14SM", "latentVectors15SM","latentVectors16SM","latentVectors17SM", "latentVectors18SM", "latentVectors19SM", "latentVectors20SM",
    #                              "latentVectors21SM","latentVectors22SM", "latentVectors23SM", "latentVectors24SM", "latentVectors25SM","latentVectors26SM","latentVectors27SM", "latentVectors28SM", "latentVectors29SM", "latentVectors30SM",
    #                              "numberOfSMNeighborAssociations", 
    #                              "numberOfMiRNANeighborAssociations", "betweennessCentralityMiRNAInMSA", "closenessCentralityMiRNAInMSA", "eigenVectorCentralityMiRNAInMSA", 
    #                              "pageRankMiRNAInMSA", "betweennessCentralitySMInMSA", "closenessCentralitySMInMSA", "eigenVectorCentralitySMInMSA",
    #                              "pageRankSMInMSA")
        colnames(loocvFeature3) <- c("latentVectors1MiRNA", "latentVectors2MiRNA", "latentVectors3MiRNA", "latentVectors4MiRNA", "latentVectors5MiRNA",
                                "latentVectors1SM","latentVectors2SM", "latentVectors3SM", "latentVectors4SM", "latentVectors5SM",
                                 "numberOfSMNeighborAssociations", 
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
  localLoocvTestingFeatureVectors2 <- BuildFeatures(localLoocvTestingIndices2)
  
  return(list("loocvTrainingFeatureVectors" = loocvTrainingFeatureVectors, 
              "globalLoocvTestingFeatureVectors" = globalLoocvTestingFeatureVectors, 
              "localLoocvTestingFeatureVectors" = localLoocvTestingFeatureVectors,
              "localLoocvTestingFeatureVectors2" = localLoocvTestingFeatureVectors2))
}

############### loops for global and local loocv #####################################################
# record the miRNA count, SM count and known MSA count
s = ALL_SM # 831
m <- max(knownMSA$V2) # 541
noOfKnownMSA <- nrow(knownMSA) # num of associations

# selection of xgboost parameters
#parameters <- list(eta = 0.2, maxDepth = 6, lambda = 1, gamma = 0,min_child_weight=3,max_delta_step=2)

# placeholder for loocv rankings
rankings <- matrix(nrow = 0, ncol = 3)
colnames(rankings) <- c("globalRankings", "localRankings_SM","localRankings_miRNA")


# loocv loops
# unit_loocv<-function(negated){
    # find negated miRNA, SM and their association's index
for(negated in 1:664){
negatedMiRNA <- knownMSA$V2[negated]
  negatedSM <- knownMSA$V1[negated]
  negatedIndex <- (negatedMiRNA - 1) * s + negatedSM
  #############################################
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
  ##########################################
  # randomly select 663 negative samples
  knownMSAIndices <- which(t(originalMSA) == 1, arr.ind = F)
  allIndices <- 1 : (m * s)
  
############
#k-means for re-sampling
###################

#construct the similarity vector for clustering
# negaSampleCand=allIndices[-knownMSAIndices]
# initialData_m=matrix(rep(0,length(negaSampleCand)*(m-1)),nrow=length(negaSampleCand),ncol=m-1)
# initialData_s=matrix(rep(0,length(negaSampleCand)*(s-1)),nrow=length(negaSampleCand),ncol=s-1)
# count=0
# for (cand in negaSampleCand){
#     i <- ifelse(cand %% ALL_SM == 0, cand / ALL_SM, as.integer(cand/ ALL_SM) + 1) 
#     j <- ifelse(cand %% ALL_SM == 0, ALL_SM, cand %% ALL_SM)
#     count=count+1
#     initialData_m[count,]=similaritiesOfMiRNA[i,-i]  
#     initialData_s[count,]=similaritiesOfSM[-j,j]   
# }
# initialData=cbind(initialData_m,initialData_s)
#choose the best number for clusters


#k-medoid is way too slow, drop it
#fviz_nbclust(dataset, kmeans, method = "wss") +geom_vline(xintercept = 3, linetype = 2)
#ClusteredNegSample<-pamk(initialData,krange=2:10, diss=FALSE)
#ClusteredNegSample<-pamk(initialData,krange=2:10,criterion="asw", usepam=FALSE,scaling=FALSE, alpha=0.001, diss=FALSE,critout=FALSE)
#ClusteredNegSample<-pam(initialData,2,diss=FALSE)
#clusplot(ClusteredNegSample$pamobject,shade = TRUE)
#clusplot(pam(initialData,ClusteredNegSample$nc))
################no resampling#####################
  set.seed(666)
  negativeSampleIndices <- sample(allIndices[-knownMSAIndices], size = 663, replace = F)
  
  # find indices for training data
  positiveAndNegativeIndices <- c(knownMSAIndices, negativeSampleIndices)
  
  # find indices for global and local testing data
  globalLoocvTestingIndices <- (1 : (m * s))[-knownMSAIndices]
  negatedIndexInGlobalTesting <- which(globalLoocvTestingIndices == negatedIndex)
  negatedIndexInLocalTesting <-  which(which(originalMSA[,negatedSM] == 0) == negatedMiRNA)
  localLoocvTestingIndices <- (which(originalMSA[,negatedSM] == 0) - 1) * s + negatedSM
  negatedIndexInLocalTesting2 <-  which(which(originalMSA[negatedMiRNA,] == 0) == negatedSM)
  localLoocvTestingIndices2 <- (which(originalMSA[negatedMiRNA,] == 0) - 1) * m + negatedMiRNA
# ##############
 ############10times Re-sampling#################
# ##############
#   negativeSampleIndices<-matrix(rep(0,length(knownMSAIndices)*10),nrow=663,ncol=10)
#   # negativeSampleIndices<-matrix(rep(0,80000),nrow=8000,ncol=10)
#   positiveAndNegativeIndices<-matrix(rep(0,length(knownMSAIndices)*20),nrow=10,ncol=1326)
#   # positiveAndNegativeIndices<-matrix(rep(0,(length(knownMSAIndices)+8000)*10),nrow=10,ncol=8663)
#   for(k in 1:10){
#   randomSeed=k*100+k*10+k
#   set.seed(randomSeed)
#   negativeSampleIndices[,k] <- as.matrix(sample(allIndices[-knownMSAIndices], size = 663, replace = F))
#   # find indices for training data
#   positiveAndNegativeIndices[k,] <- c(knownMSAIndices, negativeSampleIndices[,k])
#   }
  
#   # find indices for global and local testing data
#   globalLoocvTestingIndices <- (1 : (m * s))[-knownMSAIndices]
#   negatedIndexInGlobalTesting <- which(globalLoocvTestingIndices == negatedIndex)
#   negatedIndexInLocalTesting <- which(which(originalMSA[,negatedSM] == 0) == negatedMiRNA)
#   localLoocvTestingIndices <- (which(originalMSA[,negatedSM] == 0) - 1) * s + negatedSM         
#   negatedIndexInLocalTesting2 <- which(which(originalMSA[negatedMiRNA,] == 0) == negatedSM)
#   localLoocvTestingIndices2 <- (which(originalMSA[negatedMiRNA,] == 0) - 1) * m + negatedMiRNA  


  ########################MSA preprocessing-----fix the zeros#######################
   
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
  
  ##############################################

  # build training and testing data

  #####################10 times sampling version####################
#   for(k in 1:10){
#   trainingAndTestingData <- BuildTrainingAndTestingData(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
#                                                         negativeSampleIndices[,k], positiveAndNegativeIndices[k,], globalLoocvTestingIndices,
#                                                         localLoocvTestingIndices,localLoocvTestingIndices2)
  
#   # fit xgboost
#   xgboostLoocv <- xgboost(data = trainingAndTestingData$loocvTrainingFeatureVectors[, -1], booster = "gbtree", 
#                           label = trainingAndTestingData$loocvTrainingFeatureVectors[, 1], params = parameters, nthread = 2, nrounds = 4, 
#                           objective = "binary:logitraw")
#   if(k==1){        #declare a matrix to store all 10 weight results
#     localrklen1=nrow(trainingAndTestingData$localLoocvTestingFeatureVectors)
#     globalrklen=nrow(trainingAndTestingData$globalLoocvTestingFeatureVectors)
#     localrklen2=nrow(trainingAndTestingData$localLoocvTestingFeatureVectors2)
#     allpredictedWeightsLocal<-matrix(rep(0,localrklen1),nrow=1,ncol=localrklen1)
#     allpredictedWeightsGlobal<-matrix(rep(0,globalrklen),nrow=1,ncol=globalrklen)
#     allpredictedWeightsLocal2<-matrix(rep(0,localrklen2),nrow=1,ncol=localrklen2)
#   }
#   #get 10 different results
#     allpredictedWeightsGlobal <- rbind(allpredictedWeightsGlobal,predict(xgboostLoocv, trainingAndTestingData$globalLoocvTestingFeatureVectors))
#     allpredictedWeightsLocal <- rbind(allpredictedWeightsLocal,predict(xgboostLoocv, trainingAndTestingData$localLoocvTestingFeatureVectors))
#     allpredictedWeightsLocal2 <- rbind(allpredictedWeightsLocal2,predict(xgboostLoocv, trainingAndTestingData$localLoocvTestingFeatureVectors2))
  
#   }
# predictedWeightsGlobal<-colMeans(allpredictedWeightsGlobal)
# predictedWeightsLocal<-colMeans(allpredictedWeightsLocal)
# predictedWeightsLocal2<-colMeans(allpredictedWeightsLocal2)




####################single prediction version############################
  #prediction
 # build training and testing data
  trainingAndTestingData <- BuildTrainingAndTestingData(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
                                                        negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices,
                                                        localLoocvTestingIndices,localLoocvTestingIndices2)
  


  # # fit xgboost
  # parameters <- list(eta = 1, maxDepth = 6, lambda = 1, gamma = 0)
  # xgboostLoocv <- xgboost(data = trainingAndTestingData$loocvTrainingFeatureVectors[, -1], booster = "gbtree", 
  #                         label = trainingAndTestingData$loocvTrainingFeatureVectors[, 1], params = parameters, nthread = 2, nrounds = 2, 
  #                         objective = "binary:logitraw")

  # predictedWeightsGlobal<-predict(xgboostLoocv, trainingAndTestingData$globalLoocvTestingFeatureVectors)
  # predictedWeightsLocal<-predict(xgboostLoocv, trainingAndTestingData$localLoocvTestingFeatureVectors)
  # predictedWeightsLocal2<-predict(xgboostLoocv, trainingAndTestingData$localLoocvTestingFeatureVectors2)
  
  ################gbm##################
  X_train=trainingAndTestingData$loocvTrainingFeatureVectors[, -1]
    #Y_trian=labels
    Y_train=trainingAndTestingData$loocvTrainingFeatureVectors[, 1]


########tuned#########
#####related#####
  trControl<-trainControl(method = "cv",number = 3)
  gbmGrid <-  expand.grid(interaction.depth =4, 
                          n.trees =  400, 
                          shrinkage = c(0.06),
                          n.minobsinnode = 10)
   gbm2=caret::train(data.frame(X_train), Y_train, method = "gbm", preProcess ='scale', 
    weights = NULL,trControl=trControl,metric = "RMSE",tuneGrid = gbmGrid)
###similar###but the MiRNA not good

#####################################################using GBDT######################################
  # trControl<-trainControl(method = "cv",number = 10)
  # gbmGrid <-  expand.grid(interaction.depth =c(3), 
  #                         n.trees =  c(200,150), 
  #                         shrinkage = c(0.1),
  #                         n.minobsinnode = 10)
  #  gbm2=caret::train(data.frame(X_train), Y_train, method = "gbm", preProcess =NULL, 
  #   weights = NULL,trControl=trControl,metric = "RMSE",tuneGrid = gbmGrid)

 #       #############
 #    #use the GBDT model to predict
 #    #############
    print("Predicting scores...")
    predictedWeightsGlobal <- predict(gbm2, data.frame(trainingAndTestingData$globalLoocvTestingFeatureVectors))
    predictedWeightsLocal <- predict(gbm2, data.frame(trainingAndTestingData$localLoocvTestingFeatureVectors))
    predictedWeightsLocal2 <- predict(gbm2, data.frame(trainingAndTestingData$localLoocvTestingFeatureVectors2))

 globalRankingOfNegated <- which(sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)$ix == negatedIndexInGlobalTesting)
  localRankingOfNegated <- which(sort.int(predictedWeightsLocal, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting)
  localRankingOfNegated2 <- which(sort.int(predictedWeightsLocal2, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting2)
#####################################################USING xgboost##################################
# parameters <- list(eta = 0.02, maxDepth = 5, lambda =5, gamma =5
#           ,subsample=0.7,verbosity=2) 
#  xgboostLoocv <- xgboost(data = trainingAndTestingData$loocvTrainingFeatureVectors[,-1], booster = "gbtree", 
#                           label = trainingAndTestingData$loocvTrainingFeatureVectors[,1], params = parameters, nrounds = 100, 
#                           objective = "binary:logitraw",verbose=2)
  
#   # prediction
#   predictedWeightsGlobal <- predict(xgboostLoocv, as.matrix(trainingAndTestingData$globalLoocvTestingFeatureVectors))
#   predictedWeightsLocal <- predict(xgboostLoocv, as.matrix(trainingAndTestingData$localLoocvTestingFeatureVectors))
#  predictedWeightsLocal2 <- predict(xgboostLoocv, as.matrix(trainingAndTestingData$localLoocvTestingFeatureVectors2))

#  globalRankingOfNegated <- which(sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)$ix == negatedIndexInGlobalTesting)
#   localRankingOfNegated <- which(sort.int(predictedWeightsLocal, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting)
# localRankingOfNegated2 <- which(sort.int(predictedWeightsLocal2, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting2)


#######AUC########

 # preGlobalProb <- (predictedWeightsGlobal-min(predictedWeightsGlobal))/(max(predictedWeightsGlobal)-min(predictedWeightsGlobal))
 #  preLocalProb1<-(predictedWeightsLocal-min(predictedWeightsLocal))/(max(predictedWeightsLocal)-min(predictedWeightsLocal))

  # true_class<-str_glue("Class{PARS_human[holdout[[i]],][,151]}")
  # class_1_prob<-PARS_human_result1
  # test_set <- data.frame(obs = true_class,
  #                        Class1 = class_1_prob)
  # test_set$Class0 <- 1 - test_set$Class1
  # test_set$pred <- factor(ifelse(test_set$Class1 >= .5, "Class1", "Class0"))
  # #get the scores
  # Confu_matrix<-confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")
  # ACC<- Confu_matrix$overall["Accuracy"]
  # unitResult<-prSummary(test_set, lev = levels(test_set$obs))
  # Result[i,]<-t(as.matrix(c(unitResult,ACC)))
  ####Train AUC#####
# predictedTrain <- predict(xgboostLoocv, as.matrix(trainingAndTestingData$loocvTrainingFeatureVectors[,-1]))
predictedTrain <- predict(gbm2, as.matrix(trainingAndTestingData$loocvTrainingFeatureVectors[,-1]))
 preTrainProb <- (predictedTrain-min(predictedTrain))/(max(predictedTrain)-min(predictedTrain))
true_class<-str_glue("Class{trainingAndTestingData$loocvTrainingFeatureVectors[,1]}")
 class_1_prob<- preTrainProb
 test_set <- data.frame(obs = true_class,
                        Class1 = class_1_prob)
 test_set$Class0 <- 1 - test_set$Class1
 test_set$pred <- factor(ifelse(test_set$Class1 >= .5, "Class1", "Class0"))
 #get the scores
 Confu_matrix<-confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")
 ACC<- Confu_matrix$overall["Accuracy"]
 unitResult<-prSummary(test_set, lev = levels(test_set$obs))
 Result<-t(as.matrix(c(unitResult,ACC)))

########test Global AUC########
predictedTestG <- predict(gbm2, as.matrix(trainingAndTestingData$globalLoocvTestingFeatureVectors))
 preGlobalProb <- (predictedTestG-min(predictedTestG))/(max(predictedTestG)-min(predictedTestG))
 true_class<-str_glue("Class{0}")
 class_1_prob<- preGlobalProb
 test_set <- data.frame(obs = true_class,
                        Class1 = class_1_prob)
 test_set$Class0 <- 1 - test_set$Class1
 test_set$pred <- factor(ifelse(test_set$Class1 >= .5, "Class1", "Class0"))





  
#   X_train=trainingAndTestingData$loocvTrainingFeatureVectors[, c(-1,-15,-33,-40,-85)]
#     #Y_trian=labels
#     Y_train=trainingAndTestingData$loocvTrainingFeatureVectors[, 1]

#   trControl<-trainControl(method = "cv",number = 10)
#   gbmGrid <-  expand.grid(interaction.depth =3, 
#                           n.trees =  c(1000,1200), 
#                           shrinkage = c(0.1),
#                           n.minobsinnode = 10)
# #test maximization
#    gbm2=caret::train(data.frame(X_train), Y_train, method = "gbm", preProcess = "pca", 
#     weights = NULL,trControl=trControl,metric = "RMSE",tuneGrid = gbmGrid)

#     #############
#     #use the model to predict
#     #############
#     print("Predicting scores...")
#     predictedWeightsGlobal <- predict(gbm2, data.frame(trainingAndTestingData$globalLoocvTestingFeatureVectors[,c(-14,-32,-39,-84)]))
#     predictedWeightsLocal <- predict(gbm2, data.frame(trainingAndTestingData$localLoocvTestingFeatureVectors[,c(-14,-32,-39,-84)]))
#     predictedWeightsLocal2 <- predict(gbm2, data.frame(trainingAndTestingData$localLoocvTestingFeatureVectors2[,c(-14,-32,-39,-84)]))




  # build rankings

  ###############average_ranking###############

  # globalscoreDicrease=sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)
  # localscoreDicrease=sort.int(predictedWeightsLocal, decreasing = T, index.return = T)
  # localscoreDicrease2=sort.int(predictedWeightsLocal2, decreasing = T, index.return = T)

  # globalRankingOfNegated <- which(globalscoreDicrease$ix == negatedIndexInGlobalTesting)
  # localRankingOfNegated <- which(localscoreDicrease$ix == negatedIndexInLocalTesting)
  # localRankingOfNegated2 <- which(localscoreDicrease2$ix == negatedIndexInLocalTesting2)
  
  # globalfinalrank <- mean(which(globalscoreDicrease$x==globalscoreDicrease$x[globalRankingOfNegated]))
  # localfinalrank <- mean(which(localscoreDicrease$x==localscoreDicrease$x[localRankingOfNegated]))
  # localfinalrank2 <- mean(which(localscoreDicrease2$x==localscoreDicrease2$x[localRankingOfNegated2]))

  # rankings <- rbind(rankings, c(globalfinalrank, localfinalrank, localfinalrank2))

  #######################accurate ranking###################

 
  rankings <- rbind(rankings, c(globalRankingOfNegated, localRankingOfNegated, localRankingOfNegated2))
#   return(rankings)
# }
 
}
  



#################################use multi-cores running################################
#REMEMBER : change 3-cores to 7-cores on workstation
# clus <- makeCluster(3)

# ################################declare the function on each node###################################
# negated<-c(1:5)

# ############################################use clusters################################
# clusterExport(clus, c("MshapeOfMiRNA","ALL_MiRNA","MshapeOfSM","s","m","noOfKnownMSA","rankings",
#   "negated","knownMSA","similaritiesOfMiRNA","similaritiesOfSM","ALL_SM"),envir=environment())
# #########################clusterEvalQ_start################################
# clusterExport(clus,c("BuildTrainingAndTestingData","unit_loocv"))
# clusterEvalQ(clus,c(library(igraph),library(rNMF),library(xgboost),library(fpc),
# library(cluster),library(gbm),library(caret),library(plyr),library(parallel),library(snow)))

#   # resultt<-parLapply(clus,negated,unit_loop.function)
#   resultt<-pblapply(negated, unit_loocv,cl=clus)
#   res.df <- do.call('rbind',resultt) 
#   stopCluster(clus)
  




# write rankings to disk
#write.csv(res.df, file = "./test/old_feature_XGBOOST/tuned1-21S.csv", row.names = F)



# write rankings to disk
write.csv(rankings, file = "./test/old_feature_XGBOOST/New_xgboost.csv", row.names = F)
# write.table(rankings, file = "./test/newrgbdt.txt", col.names = T,row.names = F, sep = "\t")