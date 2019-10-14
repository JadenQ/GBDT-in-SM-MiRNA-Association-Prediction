
#################################use multi-cores running################################
#REMEMBER : change 3-cores to 7-cores on workstation
clus <- makeCluster(3)

################################declare the function on each node###################################
negated<-c(1:664)
excNegated<-c()
for(i in 1:664){
  excNegated<-rbind(excNegated,t(negated[-i]))
}
excNegated<-t(excNegated)  #663rows,664cols

clusterEvalQ(clus, unit_loop.function<-function(knownMSA,similaritiesOfMiRNA,similaritiesOfSM){
  #LOOCV loop
  s = ALL_SM # 831
  m <- max(knownMSA$V2) # 541
  noOfKnownMSA <- nrow(knownMSA) # num of associations
  rankings <- matrix(nrow = 0, ncol = 3)
  colnames(rankings) <- c("globalRankings", "localRankings_SM","localRankings_miRNA")
  



  MshapeOfMiRNA=dim(similaritiesOfMiRNA)

  ALL_MiRNA=MshapeOfMiRNA[c(1)]

  MshapeOfSM=dim(similaritiesOfSM)

  ALL_SM=MshapeOfSM[c(1)]
  # find negated miRNA, SM and their association's index
  negatedMiRNA <- knownMSA$V2[negated]
  negatedSM <- knownMSA$V1[negated]
  negatedIndex <- (negatedMiRNA - 1) * s + negatedSM
  #############################################
  # build MSA matrix
  loocvKnownMSA<-apply(excNegated,2,function(x) knownMSA[x,])
  # originalMSA <- matrix(data = rep(0, m * s), nrow = m, ncol = s)
  originalMSA<-lapply(1:664, function(x) matrix(data=rep(0,m*s), nrow=m, ncol=s))

  for(i in 1 : m) {
          print(i)
    # negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
    negatedAssociations <- lapply(loocvKnownMSA,function(x) subset(x,V2==i,select=V1))
    for(j in 1 : s) {
      # if (j %in% negatedAssociations$V1) {
      #   originalMSA[i, j] <-  1
      # }
      valueOfReac<-lapply(negatedAssociations,function(x) ifelse(j %in% x$V1,1,0))
      for(k in 1:664){
        originalMSA[[k]][i,j]<-as.numeric(valueOfReac[k])
      }
    }
  }

  MSA<-lapply(1:664, function(x) matrix(data=rep(0,m*s), nrow=s, ncol=m))
  for(i in 1 : s) {
          print(i)
    # negatedAssociations <- subset(loocvKnownMSA, V2 == i, select = V1)
    negatedAssociations <- lapply(loocvKnownMSA,function(x) subset(x,V1==i,select=V2))
    for(j in 1 : m) {
      valueOfReac<-lapply(negatedAssociations,function(x) ifelse(j %in% x$V2,1,0))
      for(k in 1:664){
       MSA[[k]][i,j]<-as.numeric(valueOfReac[k])
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
#########################
# How to release the unused variants?????????
#-negated,-knownMSAIndices,MSA!!!!!!!!!!
########################
  ###############
  #Sampling
  ###############
  # knownMSAIndices <- which(t(originalMSA) == 1, arr.ind = F)
  knownMSAIndices<-lapply(MSA, function(x) which(t(x)==1,arr.ind=F))
  
  unknownMSAIndices<-c()
  for(i in 1:length(allIndices)){
    excNegated<-rbind(excNegated,t(negated[-i]))
  }

  allIndices <- 1 : (m * s)
##########################################下一步：将下面的变量矩阵化，然后组成data.frame用apply调用下面的subGraphFeature################
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
  FeatureVOfLocaltestSamples1<-FeatureVOfGlobalTestSamples[which(globalLoocvTestingIndices%in%localLoocvTestingIndices),]
  # SampleIndices<-localLoocvTestingIndices
  # FeatureVOfLocaltestSamples1<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)
  print("Build features for local2 testing data...")
  FeatureVOfLocaltestSamples2<-FeatureVOfGlobalTestSamples[which(globalLoocvTestingIndices%in%localLoocvTestingIndices2),]
  # SampleIndices<-localLoocvTestingIndices2
  # FeatureVOfLocaltestSamples2<-subGraphFeature(MSA, similaritiesOfMiRNA, Wm,Ws,Km,Ks,KmIndex,KsIndex,similaritiesOfSM, m, s, SampleIndices)

#############
#grid search find the best model
#############
  ###############
  #Training model
  ###############
  set.seed(666)
  print("Training the model...")
  X_train=FeatureVOfTrainingSamples[, -1]
  #Y_trian=labels
  Y_train=FeatureVOfTrainingSamples[, 1]
  metric = ifelse(is.factor(Y_train), "Accuracy", "RMSE")
  trControl=trainControl(method="none")
  gbm1=caret::train(data.frame(X_train), Y_train, method = "gbm", preProcess = NULL, 
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

  return(rankings)}
)#this is for the end of clusterEvalQ
clusterExport(clus,"custom.function")



write.csv(rankings, file = "./test/1-100simi.csv", row.names = F)
write.table(rankings, file = "./test/1-100simi.txt", col.names = T,row.names = F, sep = "\t")

