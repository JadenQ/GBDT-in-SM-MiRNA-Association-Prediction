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

############
#function for feature extraction
############

# BuildTrainingAndTestingData <- function(MSA, similaritiesOfMiRNA, similaritiesOfSM, m, s, knownMSAIndices, 
#                                         negativeSampleIndices, positiveAndNegativeIndices, globalLoocvTestingIndices, localLoocvTestingIndices,localLoocvTestingIndices2) {






s = 39 # 831
m <- max(knownMSA$V2) # 541
noOfKnownMSA <- nrow(knownMSA) # num of associations
for(negated in 1 : 1) {
  # find negated miRNA, SM and their association's index
  negatedMiRNA <- knownMSA$V2[negated]
  negatedSM <- knownMSA$V1[negated]
  negatedIndex <- (negatedMiRNA - 1) * s + negatedSM
  #############################################
  # build MSA matrix
  loocvKnownMSA <- knownMSA[-negated, ]
  MSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  for(i in 1 : m) {
    negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
    for(j in 1 : s) {
      if (j %in% negatedAssociations$V1) {
        MSA[i, j] <-  1
      }
    }
  }
  
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
  
  critical=0.0001
  set.seed(666)
  Wm<-matrix(runif((ALL_MiRNA*ALL_MiRNA),min=0,max=1),ALL_MiRNA,ALL_MiRNA)
  Ws<-matrix(runif((ALL_SM*ALL_SM),min=0,max=1),ALL_SM,ALL_SM)
  Wm0<-diag(x=1,ALL_MiRNA,ALL_MiRNA)
  Distance=1
  a=0.0001                        #restart probability
  LastWm<-Wm0
  iteration=0
  while(Distance>critical){
    LastWm<-Wm
    Wm<-(1-a)*Nm*LastWm+a*Wm0
    iteration=iteration+1
    Distance=norm((Wm-LastWm),type=c("2")) #2-r norm
    print(iteration)
    print(Distance)
  }
  Ws0<-diag(x=1,ALL_SM,ALL_SM)
  Distance=1
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
      if(similaritiesOfMiRNA[i,j]>0.4){ #should consider this threshold
        #link MiRNA
        MiRNANet<-add_edges(MiRNANet,c(i,j),weight=Wm[i,j])
      }   
    }
  }
  
  for(i in 1:s){
    for(j in 1:i){
      if(similaritiesOfSM[i,j]>0.41){
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
  all<-simplify(All_net)
  
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
  
  
  ###########
  #subGraph
  ###########
  for(i in 1:m){
    for(j in 1:s){
      subGraph<-make_empty_graph(n=2,directed=FALSE)
      subGraph<-set_vertex_attr(subGraph, "name",index=1, value = paste("M",i,sep = ""),label=i)
      subGraph<-set_vertex_attr(subGraph, "name",index=2, value = paste("S",j,sep = ""),label=(m+j))
 #      subGraph<-add_vertices(subGraph,6)               #k-similar MiRNA
 # #add edge for know iteractions
 #      for(k in 1:3){                           #M1->S6,S7,S8
 #        if(MSA[i,KmIndex[i,k]]==1){
 #            subGraph<-add_edges(subGraph,c(1,(5+k))),weight=1)
 #        }
 #      }
 #      for(k in 1:3){                          #S2->M3,M4,M5
 #            if(MSA[(KsIndex[j,k]-m),j]==1){
 #             subGraph<-add_edges(subGraph,c((2+K),2),weight=1)
 #        }
 #      }
 #find m->s interactions
      countS=0
      for(k in 1:s){
        if(MSA[i,k]==1){
          countS=countS+1
          subGraph<-add_vertices(subGraph,1,name=paste("S",k,sep=""),label=(m+k))
          subGraph<-add_edges(subGraph,c(1,which(V(subGraph)$label==(m+k))),weight=1)
        }
      }
#add top3 similar SM node
      countS2=0
      for(k in 1:3){
          if(KsIndex[j,k] %in% V(subGraph)$label){
            subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==KsIndex[j,k])),weight=Ks[j,k])
          }
          else
          {
            countS2=countS2+1
            subGraph<-add_vertices(subGraph,1,name=paste("S",KsIndex[j,k],sep=""),label=KsIndex[j,k])
            subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==KsIndex[j,k])),weight=Ks[j,k])
          }
      }
  #find s->m interactions
      countM=0
      for(k in 1:m){
        if(MSA[k,j]==1){
          countM=countM+1
          subGraph<-add_vertices(subGraph,1,name=paste("M",k,sep=""),label=k)
          subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==k)),weight=1)
        }
      }
  #add top3 similar MiRNA node
      countM2=0
      for(k in 1:3){
          if(KmIndex[i,k] %in% V(subGraph)$label)
            {
              subGraph<-add_edges(subGraph,c(1,which(V(subGraph)$label==KmIndex[i,k])),weight=Km[i,k])
            }
          else{
            countM2=countM2+1
            subGraph<-add_vertices(subGraph,1,name=paste("M",KmIndex[i,k],sep=""),label=KmIndex[i,k])
            subGraph<-add_edges(subGraph,c(1,which(V(subGraph)$label==KmIndex[i,k])),weight=Km[i,k])
          }  
        }

       
  #connect all similar elements 
  ss=2+countS+countS2 #the range of SM
  mm=ss+countM+countM2 #the range of MiRNA
  #add edges for similar SM
  for(k in 3:ss){
    for(l in 3:ss)
    {
      positionI=(V(subGraph)$label[k]-m)
      positionII=(V(subGraph)$label[l]-m)
      if(similaritiesOfSM[positionI,positionII]>0.5){
        subGraph<-add_edges(subGraph,c(k,l),weight=Ws[positionI,positionII])
      }

    }
  }
  #add edges for similar MiRNAs
  for(k in (ss+1):mm){
    for(l in (ss+1):mm)
    {
      positionI=V(subGraph)$label[k]
      positionII=V(subGraph)$label[l]
      if(similaritiesOfMiRNA[positionI,positionII]>0.5){
        subGraph<-add_edges(subGraph,c(k,l),weight=Wm[positionI,positionII])
      }
    }
  }
  #add weight=1 within new added nodes
  lineS1=3+countS
  lineS2=lineS1+countS2
  lineM1=lineS2+countM+1
  lineM2=lineM1+countM2
  for(k in lineS1:lineS2){
    for (l in lineM1:lineM2) {
      if(MSA[V(subGraph)$label[l],(V(subGraph)$label[k]-m)]==1){
        subGraph<-add_edges(subGraph,c(k,l),weight=1)
      }
    }
  }
  subGraph<-simplify(subGraph)
  #find good paths
  Type1PathI<-matrix(rep(0,3),1) #the first row go with 0
  Type2PathI<-matrix(rep(0,4),1)
  all_path<-all_simple_paths(subGraph,1,2)
  for(k in 1:length(all_path)){      
      if(length(all_path[[k]])==3){
         Type1PathI<-rbind(Type1PathI,t(as.matrix(all_path[[k]])))
      }
      else if(length(all_path[[k]])==4){
        Type2PathI<-rbind(Type2PathI,t(as.matrix(all_path[[k]])))
      }
  
  #feature extraction
  # distinct 6 different paths
  C1<-matrix(rep(0,3),1)  #M,M,S
  C2<-matrix(rep(0,3),1)  #M,S,S
  C3<-matrix(rep(0,4),1)  #M,S,S,S
  C4<-matrix(rep(0,4),1)  #M,M,S,S
  C5<-matrix(rep(0,4),1)  #M,M,M,S
  C6<-matrix(rep(0,4),1)  #M,S,M,S
  for(k in 2:nrow(Type1PathI)){
    if(Type1PathI[k,2]>mm){
      C1<-rbind(C1,Type1PathI[k,])
    }
    else{
      C2<-rbind(C2,Type1PathI[k,])
    }
  }
  for(k in 2:nrow(Type2PathI)){
    if((Type1PathI[k,2]<=ss)&&(Type1PathI[k,3]<=ss)){
      C3<-rbind(C3,Type2PathI[k,])
    }
    else if((Type1PathI[k,2]>ss)&&(Type1PathI[k,3]<=ss)){
      C4<-rbind(C4,Type2PathI[k,])
    }
    else if((Type1PathI[k,2]>ss)&&(Type1PathI[k,3]>ss)){
      C5<-rbind(C5,Type2PathI[k,])
    }
    else if((Type1PathI[k,2]<=ss)&&(Type1PathI[k,3]>ss)){
      C6<-rbind(C6,Type2PathI[k,])
    }    
  }  

  #get edgeIDs
  C1edgeID<-matrix(rep(0,2*(nrow(C1)-1)),(nrow(C1)-1),2)
  C2edgeID<-matrix(rep(0,2*(nrow(C2)-1)),(nrow(C2)-1),2)
  C3edgeID<-matrix(rep(0,3*(nrow(C3)-1)),(nrow(C3)-1),3)
  C4edgeID<-matrix(rep(0,3*(nrow(C4)-1)),(nrow(C4)-1),3)
  C5edgeID<-matrix(rep(0,3*(nrow(C5)-1)),(nrow(C5)-1),3)
  C6edgeID<-matrix(rep(0,3*(nrow(C6)-1)),(nrow(C6)-1),3)

  for(k in 2:nrow(C1)){
    EP1=rep(C1[k,],each=2)[-1]
    EP1=EP1[-length(EP1)]
    C1edgeID[k,]<-get.edge.ids(subGraph,EP1)
  }
  for(k in 2:nrow(C2)){
    EP2=rep(C2[k,],each=2)[-1]
    EP2=EP2[-length(EP2)]
    C2edgeID[k,]<-get.edge.ids(subGraph,EP2)
  }
  for(k in 2:nrow(C3)){
    EP3=rep(C3[k,],each=2)[-1]
    EP3=EP3[-length(EP3)]
    C3edgeID[k,]<-get.edge.ids(subGraph,EP3)
  }
  for(k in 2:nrow(C4)){
    EP4=rep(C4[k,],each=2)[-1]
    EP4=EP4[-length(EP4)]
    C4edgeID[k,]<-get.edge.ids(subGraph,EP4)
  }
  for(k in 2:nrow(C5)){
    EP5=rep(C5[k,],each=2)[-1]
    EP5=EP5[-length(EP5)]
    C5edgeID[k,]<-get.edge.ids(subGraph,EP5)
  }
  for(k in 2:nrow(C6)){
    EP6=rep(C6[k,],each=2)[-1]
    EP6=EP6[-length(EP6)]
    C6edgeID[k,]<-get.edge.ids(subGraph,EP6)
  }

#Featrure1


  Feature1<-matrix(rep(0,6),1,6)    
  Feature2<-matrix(rep(0,6),1,6)
  Feature3<-matrix(rep(0,6),1,6)
  
      }
    }#SUBGRAPH i
  }#subgraph j
  
}
