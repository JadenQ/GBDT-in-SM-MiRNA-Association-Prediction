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

subGraphFeature<-function(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices){
#globalLoocvTestingIndices, localLoocvTestingIndices,localLoocvTestingIndices2
  

  ###########
  #subGraph
  ###########
  FVector<-matrix(rep(0,19*length(SampleIndices)),length(SampleIndices),19)
  countSamples=1
  MiRNAIndices <- ifelse(SampleIndices %% ALL_SM == 0, SampleIndices / ALL_SM, 
                                              as.integer(SampleIndices / ALL_SM) + 1) 
  SMIndices <- ifelse(SampleIndices %% ALL_SM == 0, ALL_SM, SampleIndices %% ALL_SM)
  thresholdOfSimilarity=0.8
  thresholdOfInter=0.8

for(times in 1:length(SampleIndices)){
  i<-MiRNAIndices[times]
  j<-SMIndices[times]

  subGraph<-make_empty_graph(n=2,directed=FALSE)
      subGraph<-set_vertex_attr(subGraph, "name",index=1, value = paste("M",i,sep = ""))
      subGraph<-set_vertex_attr(subGraph,"label",index=1,value=i)
      subGraph<-set_vertex_attr(subGraph, "name",index=2, value = paste("S",j,sep = ""))
      subGraph<-set_vertex_attr(subGraph,"label",index=2,value=(m+j))

      countS=0
      for(k in 1:s){
        if((MSA[i,k]>thresholdOfInter)&&(k!=j)){
          countS=countS+1
          subGraph<-add_vertices(subGraph,1,name=paste("S",k,sep=""),label=(m+k))
          subGraph<-add_edges(subGraph,c(1,which(V(subGraph)$label==(m+k))),weight=MSA[i,k])
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
            subGraph<-add_vertices(subGraph,1,name=paste("S",(KsIndex[j,k]-m),sep=""),label=KsIndex[j,k])
            subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==KsIndex[j,k])),weight=Ks[j,k])
          }
      }
  #find s->m interactions
      countM=0
      for(k in 1:m){
        if((MSA[k,j]>thresholdOfInter)&&(k!=i)){
          countM=countM+1
          subGraph<-add_vertices(subGraph,1,name=paste("M",k,sep=""),label=k)
          subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==k)),weight=MSA[k,j])
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
      if(similaritiesOfSM[positionI,positionII]>thresholdOfSimilarity){
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
      if(similaritiesOfMiRNA[positionI,positionII]>thresholdOfSimilarity){
        subGraph<-add_edges(subGraph,c(k,l),weight=Wm[positionI,positionII])
      }
    }
  }
  #add weight=1 within new added nodes
  lineS1=2+countS
  lineS2=lineS1+countS2
  lineM1=lineS2+countM
  lineM2=lineM1+countM2
  showTheSub<-V(subGraph)$label
  if((countS2!=0)&&(countM2!=0)){
    for(k in (lineS1+1):lineS2){
      for (l in (lineM1+1):lineM2) {
        positionI<-V(subGraph)$label[l]
        positionII<-V(subGraph)$label[k]-m
        if(MSA[positionI,positionII]>thresholdOfInter){
          subGraph<-add_edges(subGraph,c(k,l),weight=MSA[positionI,positionII])
        }
      }
    }
  }

  subGraph<-simplify(subGraph)
  print("Get subGraph...")
  #find good paths
  Type1PathI<-matrix(rep(0,3),1) #the first row go with 0
  Type2PathI<-matrix(rep(0,4),1)
  all_path<-all_simple_paths(subGraph,1,2)
  print("Find all paths...")
  print(countSamples)
  #smoothing
  if(length(all_path)<1){
    subGraph<-add.edges(subGraph,c(1,2),weight=0.000001)
    all_path<-all_simple_paths(subGraph,1,2)
  }
  for(k in 1:length(all_path)){      
      if(length(all_path[[k]])==3){
         Type1PathI<-rbind(Type1PathI,t(as.matrix(all_path[[k]])))
      }
      else if(length(all_path[[k]])==4){
        Type2PathI<-rbind(Type2PathI,t(as.matrix(all_path[[k]])))
      }
  }
  #feature extraction
  # distinct 6 different paths
  C1<-matrix(rep(0,3),1)  #M,M,S
  C2<-matrix(rep(0,3),1)  #M,S,S
  C3<-matrix(rep(0,4),1)  #M,S,S,S
  C4<-matrix(rep(0,4),1)  #M,M,S,S
  C5<-matrix(rep(0,4),1)  #M,M,M,S
  C6<-matrix(rep(0,4),1)  #M,S,M,S
  if(nrow(Type1PathI)>1){
      for(k in 2:nrow(Type1PathI)){
        if(Type1PathI[k,2]<mm){
          C1<-rbind(C1,Type1PathI[k,])
        }
        else{
          C2<-rbind(C2,Type1PathI[k,])
        }
    }
  }
  else{
    C1Weight<-matrix(rep(0,2),1)
    C2Weight<-matrix(rep(0,2),1)
  }
  if(nrow(Type2PathI)>1){
    for(k in 2:nrow(Type2PathI)){
      if((Type2PathI[k,2]<=ss)&&(Type2PathI[k,3]<=ss)){
        C3<-rbind(C3,Type2PathI[k,])
      }
      else if((Type2PathI[k,2]>ss)&&(Type2PathI[k,3]<=ss)){
        C4<-rbind(C4,Type2PathI[k,])
      }
      else if((Type2PathI[k,2]>ss)&&(Type2PathI[k,3]>ss)){
        C5<-rbind(C5,Type2PathI[k,])
      }
      else if((Type2PathI[k,2]<=ss)&&(Type2PathI[k,3]>ss)){
        C6<-rbind(C6,Type2PathI[k,])
      }    
    }
  }
  else{
    C3Weight<-matrix(rep(0,3),1)
    C4Weight<-matrix(rep(0,3),1)
    C5Weight<-matrix(rep(0,3),1)
    C6Weight<-matrix(rep(0,3),1)
  }
    

  if(nrow(C1)>1){
    C1edgeID<-matrix(rep(0,2*(nrow(C1)-1)),(nrow(C1)-1),2)
    for(k in 2:nrow(C1)){
        EP1=rep(C1[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C1edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
    if(nrow(C1edgeID)>1){
        C1Weight<-as.matrix(E(subGraph)$weight[C1edgeID[1,]])
      for(k in 2:nrow(C1edgeID)){
        C1Weight<-rbind(C1Weight,as.matrix(E(subGraph)$weight[C1edgeID[k,]]))
      }
    }
    else{
        C1Weight<-as.matrix(E(subGraph)$weight[C1edgeID[1,]])
    }

  }
  else {C1Weight<-matrix(rep(0,2),1)}

  if(nrow(C2)>1){
    C2edgeID<-matrix(rep(0,2*(nrow(C2)-1)),(nrow(C2)-1),2)
    for(k in 2:nrow(C2)){
        EP1=rep(C2[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C2edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
    if(nrow(C2edgeID)>1){
        C2Weight<-as.matrix(E(subGraph)$weight[C2edgeID[1,]])
      for(k in 2:nrow(C2edgeID)){
        C2Weight<-rbind(C2Weight,as.matrix(E(subGraph)$weight[C2edgeID[k,]]))
      }
    }
    else{
        C2Weight<-as.matrix(E(subGraph)$weight[C2edgeID[1,]])
    }

  }
  else {C2Weight<-matrix(rep(0,2),1)}


  if(nrow(C3)>1){
    C3edgeID<-matrix(rep(0,3*(nrow(C3)-1)),(nrow(C3)-1),3)
    for(k in 2:nrow(C3)){
        EP1=rep(C3[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C3edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
      if(nrow(C3edgeID)>1){
        C3Weight<-as.matrix(E(subGraph)$weight[C3edgeID[1,]])
      for(k in 2:nrow(C3edgeID)){
        C3Weight<-rbind(C3Weight,as.matrix(E(subGraph)$weight[C3edgeID[k,]]))
      }
    }
    else{
        C3Weight<-as.matrix(E(subGraph)$weight[C3edgeID[1,]])
    }

  }
  else {C3Weight<-matrix(rep(0,3),1)}

  if(nrow(C4)>1){
    C4edgeID<-matrix(rep(0,3*(nrow(C4)-1)),(nrow(C4)-1),3)
    for(k in 2:nrow(C4)){
        EP1=rep(C4[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C4edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
    if(nrow(C4edgeID)>1){
        C4Weight<-as.matrix(E(subGraph)$weight[C4edgeID[1,]])
      for(k in 2:nrow(C4edgeID)){
        C4Weight<-rbind(C4Weight,as.matrix(E(subGraph)$weight[C4edgeID[k,]]))
      }
    }
    else{
        C4Weight<-as.matrix(E(subGraph)$weight[C4edgeID[1,]])
    }
      

  }
  else {C4Weight<-matrix(rep(0,3),1)} 
  #C5
  if(nrow(C5)>1){
    C5edgeID<-matrix(rep(0,3*(nrow(C5)-1)),(nrow(C5)-1),3)
    for(k in 2:nrow(C5)){
        EP1=rep(C5[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C5edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
    if(nrow(C5edgeID)>1){
        C5Weight<-as.matrix(E(subGraph)$weight[C5edgeID[1,]])
      for(k in 2:nrow(C5edgeID)){
        C5Weight<-rbind(C5Weight,as.matrix(E(subGraph)$weight[C5edgeID[k,]]))
      }
    }
    else{
        C5Weight<-as.matrix(E(subGraph)$weight[C5edgeID[1,]])
    }

  }
  else {C5Weight<-matrix(rep(0,3),1)} 
#C6
  if(nrow(C6)>1){
    C6edgeID<-matrix(rep(0,3*(nrow(C6)-1)),(nrow(C6)-1),3)
    for(k in 2:nrow(C6)){
        EP1=rep(C6[k,],each=2)[-1]
        EP1=EP1[-length(EP1)]
        C6edgeID[(k-1),]<-get.edge.ids(subGraph,EP1)
      }
    #format the weights
    if(nrow(C6edgeID)>1){
        C6Weight<-as.matrix(E(subGraph)$weight[C6edgeID[1,]])
      for(k in 2:nrow(C6edgeID)){
        C6Weight<-rbind(C6Weight,as.matrix(E(subGraph)$weight[C6edgeID[k,]]))
      }
    }
    else{
        C6Weight<-as.matrix(E(subGraph)$weight[C6edgeID[1,]])
    }

  }

  else {C6Weight<-matrix(rep(0,3),1)}




#Featrures
  Feature1<-matrix(rep(0,6),1,6) 
  Feature2<-matrix(rep(0,6),1,6)
  Feature3<-matrix(rep(0,6),1,6)   

  if(nrow(C1Weight)>1){
    pro<-as.matrix(prod(C1Weight[1,]))
    for(k in 2:nrow(C1Weight)){
      pro<-rbind(pro,prod(C1Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C1Weight[1,]))
  }

  Feature1[1]<-colSums(pro)
  Feature2[1]<-max(pro)
  Feature3[1]<-(nrow(C1)-1)

  if(nrow(C2Weight)>1){
    pro<-as.matrix(prod(C2Weight[1,]))
    for(k in 2:nrow(C2Weight)){
      pro<-rbind(pro,prod(C2Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C2Weight[1,]))
  }
  Feature1[2]<-colSums(pro)
  Feature2[2]<-max(pro)
  Feature3[2]<-(nrow(C2)-1)

  if(nrow(C3Weight)>1){
    pro<-as.matrix(prod(C3Weight[1,]))
    for(k in 2:nrow(C3Weight)){
      pro<-rbind(pro,prod(C3Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C3Weight[1,]))
  }
  Feature1[3]<-colSums(pro)
  Feature2[3]<-max(pro)
  Feature3[3]<-(nrow(C3)-1)

  if(nrow(C4Weight)>1){
    pro<-as.matrix(prod(C4Weight[1,]))
    for(k in 2:nrow(C4Weight)){
      pro<-rbind(pro,prod(C4Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C4Weight[1,]))
  }
  Feature1[4]<-colSums(pro)
  Feature2[4]<-max(pro)
  Feature3[4]<-(nrow(C4)-1)

  if(nrow(C5Weight)>1){
    pro<-as.matrix(prod(C5Weight[1,]))
    for(k in 2:nrow(C5Weight)){
      pro<-rbind(pro,prod(C5Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C5Weight[1,]))
  }
  Feature1[5]<-colSums(pro)
  Feature2[5]<-max(pro)
  Feature3[5]<-(nrow(C5)-1)

   if(nrow(C6Weight)>1){
    pro<-as.matrix(prod(C6Weight[1,]))
    for(k in 2:nrow(C6Weight)){
      pro<-rbind(pro,prod(C6Weight[k,]))
   }
  }
  else{
    pro<-as.matrix(prod(C6Weight[1,]))
  }
  Feature1[6]<-colSums(pro)
  Feature2[6]<-max(pro)
  Feature3[6]<-(nrow(C6)-1)

  if(originalMSA[i,j]==1){
    Label=1
  }
  else{Label=0}
  Label<-as.matrix(Label)

  FVector[countSamples,]<-cbind(Label,Feature1,Feature2,Feature3)
  print("Get feature Vectors...")
  print(countSamples)
  countSamples=countSamples+1
}# End the subgraph 
  print("Get feature Vectors Matrix...")
  return(FVector)       #return a matrix composed of Feature Vectors from Training or Testing samples
}

#LOOCV loop
s = ALL_SM # 831
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
  originalMSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  MSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  for(i in 1 : m) {
    negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
    for(j in 1 : s) {
      if (j %in% negatedAssociations$V1) {
        originalMSA[i, j] <-  1
      }
    }
  }
  MSA <- originalMSA
#FIX the zero
 for(i in 1:m){
    for(j in 1:s){
      if(originalMSA[i,j]==0){
        ##find similar MiRNA by LINEs
        SimilarityDeacreaseOfMiRNA=sort.int(similaritiesOfMiRNA[i,],decreasing=T,index.return=T)
        count1=0                    #count the top similar MiRNAs
        count2=0
        similarMiRNAIndex<-rep(0,3) #preset index of top3 similar miRNA
        similarMiRNASValue<-rep(0,3)#preset similarity value of top3 similar miRNA
        for(k in 1:m){                             #find similar MiRNA by LINEs
            # flag=ifelse(SimilarityDeacreaseOfMiRNA$x[k]>0.5,1,-1)
            # while(flag==1&&count1<3){               #only find top3 similar MiRNAs with similarity value above 0.5
              if(originalMSA[k,j]==1&&count1<3){      #find top3 similar SM in all other MiRNAs
               count1=count1+1
               similarMiRNAIndex[count1]=SimilarityDeacreaseOfMiRNA$ix[k]
               similarMiRNASValue[count1]=SimilarityDeacreaseOfMiRNA$x[k]
              }
              else next
          }
          predictSVsepMiRNA<-rep(0,3)
          for(l in 1:3){
            if(is.nan(similarMiRNASValue[l])){similarMiRNASValue[l]=0}
            predictSVsepMiRNA[l]=(similarMiRNASValue[l]^2)/sum(similarMiRNASValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
            if(is.nan(predictSVsepMiRNA[l])){predictSVsepMiRNA[l]=0}
            else next
          }
          predictSVMiRNA=sum(predictSVsepMiRNA)
          ##find similar SM by ROWs##
          SimilarityDeacreaseOfSM=sort.int(similaritiesOfSM[,j],decreasing=T,index.return=T)
          similarSMIndex<-rep(0,3) #preset index of top3 similar SM
          similarSMSValue<-rep(0,3)#preset similarity value of top3 similar SM
          for(k in 1:s){                           #find similar SM by ROWs
           # flag=ifelse(SimilarityDeacreaseOfSM$x[k]>0.5,1,-1)
           # while(flag==1&&count2<3){               #only find top3 similar SM with similarity value above 0.5                                
              if(originalMSA[i,k]==1&&count2<3){      #find top3 similar SM in all other SMs
               count2=count2+1
               similarSMIndex[count2]=SimilarityDeacreaseOfSM$ix[k]
               similarSMSValue[count2]=SimilarityDeacreaseOfSM$x[k]
              }
              else next
            }
    #}    
          predictSVsepSM<-rep(0,3)
          for(l in 1:3){
            if(is.nan(similarSMSValue[l])){similarSMSValue[l]=0}
            predictSVsepSM[l]=(similarSMSValue[l]^2)/sum(similarSMSValue) # predict simi value=(sv1)*w1+(sv2)*w2+(sv3)*w3  P.S. wi=svi/sum[sv]
            if(is.nan(predictSVsepSM[l])){predictSVsepSM[l]=0}
            else next
          }
          predictSVSM=sum(predictSVsepSM)
          MSA[i,j]=(predictSVMiRNA+predictSVSM)/2
        }
        
        else next
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
    KmIndex[i,]<-t(as.matrix(sort.int(Wm[i,],decreasing = T,index.return = T)$ix[2:4]))
  }
  for(i in 1:s){
    Ks[i,]<-t(as.matrix(sort.int(Ws[i,],decreasing = T,index.return = T)$x[2:4]))
    KsIndex[i,]<-addForGraph+t(as.matrix(sort.int(Ws[i,],decreasing = T,index.return = T)$ix[2:4]))
  }

  ###############
  #Sampling
  ###############
  knownMSAIndices <- which(t(originalMSA) == 1, arr.ind = F)
  allIndices <- 1 : (m * s)

  negativeSampleIndices<-sample(allIndices[-knownMSAIndices], size = 663, replace = F)
  positiveAndNegativeIndices <- c(knownMSAIndices, negativeSampleIndices)
  #find global and local Indexs
  globalLoocvTestingIndices <- (1 : (m * s))[-knownMSAIndices]
  negatedIndexInGlobalTesting <- which(globalLoocvTestingIndices == negatedIndex)
  negatedIndexInLocalTesting <- which(which(originalMSA[,negatedSM] == 0) == negatedMiRNA)
  localLoocvTestingIndices <- (which(originalMSA[,negatedSM] == 0) - 1) * s + negatedSM         
  negatedIndexInLocalTesting2 <- which(which(originalMSA[negatedMiRNA,] == 0) == negatedSM)
  localLoocvTestingIndices2 <- (which(originalMSA[negatedMiRNA,] == 0) - 1) * m + negatedMiRNA 

  ################
  #Use the Feature extraction Function
  ################
  print("Build features for the training data...")
  SampleIndices<-positiveAndNegativeIndices
  FeatureVOfTrainingSamples<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)
  print("Build features for global testing data...")
  SampleIndices<-globalLoocvTestingIndices
  FeatureVOfGlobalTestSamples<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)
  print("Build features for local1 testing data...")
  SampleIndices<-localLoocvTestingIndices
  FeatureVOfLocaltestSamples1<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)
  print("Build features for local2 testing data...")
  SampleIndices<-localLoocvTestingIndices2
  FeatureVOfLocaltestSamples2<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)


  ###############
  #Training model
  ###############
  print("Training the model...")
  X_train=FeatureVOfTrainingSamples[, -1]
  #Y_trian=labels
  Y_train=FeatureVOfTrainingSamples[, 1]
  metric = ifelse(is.factor(Y_train), "Accuracy", "RMSE")
  trControl=trainControl(method="none")
  gbm1=train(X_train, Y_train, method = "gbm", preProcess = NULL, 
  weights = NULL,trControl=trainControl(method="none") ,metric = ifelse(is.factor(Y_train), "Accuracy", "RMSE"),
  maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE,
  TRUE),tuneGrid = NULL,tuneLength = ifelse(trControl$method == "none", 1, 3))

  #############
  #use the model to predict
  #############
  print("Predicting scores...")
  predictedWeightsGlobal <- predict(gbm1, data.frame(FeatureVOfGlobalTestSamples[,-1]))
  predictedWeightsLocal <- predict(gbm1, data.frame(FeatureVOfLocaltestSamples1[,-1]))
  predictedWeightsLocal2 <- predict(gbm1, data.frame(FeatureVOfLocaltestSamples2[,-1]))
  
  ##########
  #get the ranking
  ##########
  print("Give the ranking...")
  globalRankingOfNegated <- which(sort.int(predictedWeightsGlobal, decreasing = T, index.return = T)$ix == negatedIndexInGlobalTesting)
  localRankingOfNegated <- which(sort.int(predictedWeightsLocal, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting)
  localRankingOfNegated2 <- which(sort.int(predictedWeightsLocal2, decreasing = T, index.return = T)$ix == negatedIndexInLocalTesting2)
  rankings <- rbind(rankings, c(globalRankingOfNegated, localRankingOfNegated, localRankingOfNegated2))



  }#negated leave one

write.csv(rankings, file = "./test/newrgbdt-.csv", row.names = F)
write.table(rankings, file = "./test/newrgbdt.txt", col.names = T,row.names = F, sep = "\t")