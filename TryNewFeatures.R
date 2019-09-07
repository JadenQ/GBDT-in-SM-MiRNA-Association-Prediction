library(igraph)
library(rNMF)
library(xgboost)
library(fpc)
library(cluster)
library(gbm)
library(caret)
library(plyr)
# # related_group
# # read known miRNA-SM association dataset
knownMSA <- read.table(file = "./SM-miRNA/related/SM-miRNA_Num_A.csv", header = F,sep=",")
# read miRNA functional similarity matrix
similaritiesOfMiRNA <- as.matrix(read.table(file = "./SM-miRNA/related/miRNA_similarity_matrix2.csv", header = F,sep=","))
# read SM similarity matrix
similaritiesOfSM <- as.matrix(read.table(file = "./SM-miRNA/related/SM_similarity_matrix2.csv", header = F,sep=","))

MshapeOfMiRNA=dim(similaritiesOfMiRNA)

ALL_MiRNA=MshapeOfMiRNA[c(1)]

MshapeOfSM=dim(similaritiesOfSM)

ALL_SM=MshapeOfSM[c(1)]

m=286
s=39
############
#function for feature extraction
############

# BuildTrainingAndTestingData <- function(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
#                                         negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices, localLoocvTestingIndices,localLoocvTestingIndices2) {

  ##########################
  #random walk
  ##########################
  #transition probability matrix:Nm&Ns
  Nm<-matrix(rep(0,ALL_MiRNA*ALL_MiRNA),nrow=ALL_MiRNA,ncol=ALL_MiRNA)
  Ns<-matrix(rep(0,ALL_SM*ALL_SM),nrow=ALL_SM,ncol=ALL_SM)
  MiRNArowsum=rowSums(similaritiesOfMiRNA)
  SMrowsum=rowSums(similaritiesOfSM)
  for(i in 1:m){
    Nm[i,]<-similaritiesOfMiRNA[i,]/MiRNArowsum
  }
  for(i in 1:s){
    Ns[i,]<-similaritiesOfSM[i,]/SMrowsum
  }
  #RW probability matrix

  critical=0.000001
  set.seed(666)
  Wm<-matrix(runif((ALL_MiRNA*ALL_MiRNA),min=0,max=1),ALL_MiRNA,ALL_MiRNA)
  Ws<-matrix(runif((ALL_SM*ALL_SM),min=0,max=1),ALL_SM,ALL_SM)
  Wm0<-diag(x=1,ALL_MiRNA,ALL_MiRNA)
  Distance=1
  a=0.05                                #restart probability
  LastWm<-Wm0
  while(Distance>critical){
    LastWm<-Wm
    Wm<-(1-a)*Nm*LastWm+a*Wm0
    Distance=norm((Wm-LastWm),type=c("2")) #2-r norm
  }
  Ws0<-diag(x=1,ALL_SM,ALL_SM)
  Distance=1                                #renew the distance
  LastWs<-Ws0
  while(Distance>critical){
    LastWs<-Ws
    Ws<-(1-a)*Ns*LastWs+a*Ws0
    Distance=norm((Ws-LastWs),type=c("2")) #2-r norm
  }

  ###################################
  #create the network with new weight
  ###################################

  #similarity network
  Mirna_net <- make_empty_graph(n=ALL_MiRNA,directed = FALSE)
  MiRNANet <-set.vertex.attribute(Mirna_net,"name",value = paste("M",1:ALL_MiRNA,sep = "")) 
  Sm_net <- make_empty_graph(n=ALL_SM,directed=FALSE)
  SMNet <- set.vertex.attribute(Sm_net,"name",value = paste("S",1:ALL_SM,sep=""))

  for(i in 1:m){
    for(j in 1:i){
      if(similaritiesOfMiRNA[i,j]>0.5){
        #link MiRNA
        MiRNANet<-add_edges(MiRNANet,c(i,j),weight=Wm[i,j])
      }   
    }
  }

  for(i in 1:s){
    for(j in 1:i){
      if(similaritiesOfSM[i,j]>0.5){
        #link SMs
        SMNet<-add_edges(SMNet,c(i,j),weight=Ws[i,j])
      }
    }
  }
  #286+39,541+831
  All_net<-MiRNANet+SMNet
  #interaction network
  for(i in 1:m){
    for(j in 1:s){
      if(MSA[i,j]==1){
        All_net<-add_edges(All_net,c(i,m+j),weight=1)
      }
    }
  }
  All_net<-simplify(All_net)
  #########################
  #create k-similar matrix
  #########################
  Km<-matrix(rep(0,(ALL_MiRNA*3)),nrow=ALL_MiRNA,ncol=3)
  Ks<-matrix(rep(0,(ALL_SM*3)),nrow=ALL_SM,ncol=3)
  KmIndex<-matrix(rep(0,(ALL_MiRNA*3)),nrow=ALL_MiRNA,ncol=3) #index in the graph
  KsIndex<-matrix(rep(0,(ALL_SM*3)),nrow=ALL_SM,ncol=3)
  addForGraph<-matrix(rep(ALL_MiRNA,3),1,3)
  for(i in 1:m){
    Km[i,]<-t(as.matrix(sort.int(Wm[i,],decreasing = T,index.return = T)$x[2:4])) # not start at 1 because of the simi of their own
    KmIndex<-t(as.matrix(sort.int(Wm[i,],decreasing = T,index.return = T)$ix[2:4]))
  }
  for(i in 1:s){
    Ks[i,]<-t(as.matrix(sort.int(Ws[i,],decreasing = T,index.return = T)$x[2:4]))
    KsIndex<-addForGraph+t(as.matrix(sort.int(Ws[i,],decreasing = T,index.return = T)$ix[2:4]))
  }

  ################
  #micro network
  ################

  #################
  #micro network
  #too much computation
  #################
  # for(i in 1:m){
  #   for(j in 1:s){
  #     all_path<-all_simple_paths(All_net,i,(m+j))
  #     Type1PathI<-matrix(rep(0,3),1) #the first row go with 0
  #     Type2PathI<-matrix(rep(0,4),1)
  #     #all smooth: get all paths with length of 3,4 for this interaction
  #     for(k in 1:length(all_path)){      
  #       if(length(all_path[[k]])==3){
  #         Type1PathI<-rbind(Type1PathI,t(as.matrix(all_path[[k]])))
  #       }
  #       else if(length(all_path[[k]])==4){
  #         Type2PathI<-rbind(Type2PathI,t(as.matrix(all_path[[k]])))
  #       }
  #     }
  #     #find required similarity path 
  #     for(k in 1:nrow(Type1PathI)){
  #       if(Type1PathI[k,1]){
  #         #find the compatible elements
  #       }
  #     }
      
  #   }
  # }


}
#loocv loop
####################################
 
  # s = 39 # 831
  # m <- max(knownMSA$V2) # 541
  # noOfKnownMSA <- nrow(knownMSA) # num of associations
  # for(negated in 1 : 2) {
  #   # find negated miRNA, SM and their association's index
  #   negatedMiRNA <- knownMSA$V2[negated]
  #   negatedSM <- knownMSA$V1[negated]
  #   negatedIndex <- (negatedMiRNA - 1) * s + negatedSM
  #   #############################################
  #   # build MSA matrix
  #   loocvKnownMSA <- knownMSA[-negated, ]
  #   MSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  #   for(i in 1 : m) {
  #     negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
  #     for(j in 1 : s) {
  #       if (j %in% negatedAssociations$V1) {
  #         MSA[i, j] <-  1
  #       }
  #     }
  #   }
  #   miniMSA<-MSA[1:6,1:6]
  #   MSAGraph <- graph_from_incidence_matrix(incidence = miniMSA, directed = F, mode = "total")
  #   MSAGraph2<-graph_from_adjacency_matrix(miniMSA,weighted =NULL ,mode="undirected")
  #   miniMiRNA<-similaritiesOfMiRNA[1:10,1:10]
  #   miniSM<-similaritiesOfSM[1:10,1:10]
  #   MsimiGraph<-graph_from_adjacency_matrix(miniMiRNA,weighted = TRUE,
                                           
  #                                           mode="undirected")
  # }
  
  # plot(MsimiGraph, vertex.color=rainbow(10),vertex.size=30,edge.width=2,edge.color="black")
  # # build MSA graph
  # MSAGraph <- graph_from_incidence_matrix(incidence = miniMSA, directed = F, mode = "total")

  # betweennessCentralityMSA <- betweenness(MSAGraph, directed = F, normalized = T)
  # betweennessCentralityMiRNAInMSA <- betweennessCentralityMSA[1:ALL_MiRNA]
  # betweennessCentralitySMInMSA <- betweennessCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  # closenessCentralityMSA <- closeness(MSAGraph, mode = "all")
  # closenessCentralityMiRNAInMSA <- closenessCentralityMSA[1:ALL_MiRNA]
  # closenessCentralitySMInMSA <- closenessCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  # eigenVectorCentralityMSA <- eigen_centrality(MSAGraph, directed = F)$vector
  # eigenVectorCentralityMiRNAInMSA <- eigenVectorCentralityMSA[1:ALL_MiRNA]
  # eigenVectorCentralitySMInMSA <- eigenVectorCentralityMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]
  # pageRankMSA <- page.rank(MSAGraph, directed = F)$vector
  # pageRankMiRNAInMSA <- pageRankMSA[1:ALL_MiRNA]
  # pageRankSMInMSA <- pageRankMSA[(ALL_MiRNA+1):(ALL_SM+ALL_MiRNA)]

   # # number of associations between an miRNA and a SM's neighbors
  # numberOfSMNeighborAssociations <- c(rep(0, m * s))
  # for(i in 1 : m) {
  #   for(j in 1 : s) {
  #     numberOfAssociations = ifelse(MSA[i, j] == 1, -1, 0)
  #     SMNeighbors = which(t(similaritiesOfSM[j, ]) >= meanSimilaritySM, arr.ind = F)
  #     for(k in 1 : length(SMNeighbors)) {
  #       if(MSA[i, SMNeighbors[k]] == 1) {
  #         numberOfAssociations = numberOfAssociations + 1
  #       }
  #     }
  #     numberOfSMNeighborAssociations[(i-1)*s+j] <- numberOfAssociations
  #   }
  # }
  
  # # number of associations between a SM and an miRNA's neighbors
  # numberOfMiRNANeighborAssociations <- c(rep(0, m * s))
  # for(i in 1 : s) {
  #   for(j in 1 : m) {
  #     numberOfAssociations = ifelse(MSA[j, i] == 1, -1, 0)
  #     miRNANeighbors = which(t(similaritiesOfMiRNA[j, ]) >= meanSimilarityMiRNA, arr.ind = F)
  #     for(k in 1 : length(miRNANeighbors)) {
  #       if(MSA[miRNANeighbors[k]] == 1) {
  #         numberOfAssociations = numberOfAssociations + 1
  #       }
  #     }
  #     numberOfMiRNANeighborAssociations[(i-1)*m+j] <- numberOfAssociations
  #   }
  # }
  
