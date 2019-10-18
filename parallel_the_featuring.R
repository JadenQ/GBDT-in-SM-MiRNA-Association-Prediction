
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
        if((originalMSA[i,k]>thresholdOfInter)&&(k!=j)){
          countS=countS+1
          subGraph<-add_vertices(subGraph,1,name=paste("S",k,sep=""),label=(m+k))
          subGraph<-add_edges(subGraph,c(1,which(V(subGraph)$label==(m+k))),weight=originalMSA[i,k])
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
        if((originalMSA[k,j]>thresholdOfInter)&&(k!=i)){
          countM=countM+1
          subGraph<-add_vertices(subGraph,1,name=paste("M",k,sep=""),label=k)
          subGraph<-add_edges(subGraph,c(2,which(V(subGraph)$label==k)),weight=originalMSA[k,j])
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
        if(originalMSA[positionI,positionII]>thresholdOfInter){
          subGraph<-add_edges(subGraph,c(k,l),weight=originalMSA[positionI,positionII])
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
  print(times)
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

  FVector[times,]<-cbind(Label,Feature1,Feature2,Feature3)
  print("Get feature Vectors...")
  print(negated)
  print("***********")
  print(times)
}# End the subgraph 
  print("Get feature Vectors Matrix...")
  return(FVector)       #return a matrix composed of Feature Vectors from Training or Testing samples
}