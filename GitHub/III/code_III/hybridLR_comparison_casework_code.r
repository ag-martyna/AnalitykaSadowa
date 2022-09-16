source("../BW_function.R")
source("../LR_Nor_comparison_function.R")
source("../LR_KDE_comparison_function.R")
source("../rMANOVA_general_function.R")

if (!require("mvtnorm")) install.packages("mvtnorm")

population = read.table("../spectra_database.txt", header = TRUE)
data.analysed = population
data.recovered = read.table("../spectra_recovered.txt", header = TRUE) 
data.control = read.table("../spectra_control.txt", header = TRUE) 

M = m = length(unique(data.analysed$Name))
n.all=length(unique(data.analysed$Piece))
m.control = length(unique(data.control$Name))
m.recovered = length(unique(data.recovered$Name))
n.control=length(unique(data.control$Piece))
n.recovered=length(unique(data.recovered$Piece))

validationSetsNo = 5 

output.matrix.Nor = output.matrix.KDE = array(NA,c(m.recovered,m.control,validationSetsNo),dimnames=list(as.character(unique(data.recovered$Name)),as.character(unique(data.control$Name)),1:validationSetsNo))

for (s in 1:validationSetsNo)
{
  tvec = NULL
  for (g in 1:M)
  {tvec = c(tvec,sample(rep(c(1,2),each=n.all/2)))}
  
  A.1 = data.analysed[which(tvec == 1),] ##etap 1
  A.1[,"Piece"] = rep(1:(n.all/2),times=M)
  A.2 = data.analysed[which(tvec == 2),] ##etap 1
  A.2[,"Piece"] = rep(1:(n.all/2),times=M)
  
  tvec = NULL
  for (g in 1:m.control)
  {tvec = c(tvec,sample(rep(c(1,2),each=n.control/2)))}
  
  y.2.1 = data.control[which(tvec == 1),] ##etap 2
  y.2.1[,"Piece"] = rep(1:(n.control/2),times=m.control)
  y.2.2 = data.control[which(tvec == 2),] ##etap 2
  y.2.2[,"Piece"] = rep(1:(n.control/2),times=m.control)
  
  tvec = NULL
  for (g in 1:m.recovered)
  {tvec = c(tvec,sample(rep(c(1,2),each=n.recovered/2)))}
  
  y.1.1 = data.recovered[which(tvec == 1),] ##etap 2
  y.1.1[,"Piece"] = rep(1:(n.recovered/2),times=m.recovered)
  y.1.2 = data.recovered[which(tvec == 2),] ##etap 2
  y.1.2[,"Piece"] = rep(1:(n.recovered/2),times=m.recovered)
  
  for (i in 1:m.control)
  {  
    B = rbind(A.1,y.2.1[which(y.2.1$Item==i),]) ##etap 3
    B[,"Name"] = c(A.1[,3],rep("control",times=n.control/2))
    B[,"Item"] = c(rep(1:(M),each=n.all/2),rep(M+1,times=n.control/2))
    B[,"Piece"] = c(rep(1:(n.all/2),times=M),rep(1:(n.control/2),times=1))
    
    rMANOVA_res = rMANOVA_general_function(data=as.data.frame(B)) ##etap 4
    eigen_vectors = rMANOVA_res$results
    center = rMANOVA_res$center
    
    for (j in 1:m.recovered) 
    {  
      A = rbind(A.2,y.2.2[which(y.2.2$Item==i),],y.1.2[which(y.1.2$Item==j),]) ##etap 5
      A[,"Name"] = c(A.2[,1],rep("control",times=n.control/2),rep("recovered",times=n.recovered/2))
      A[,"Item"] = c(rep(1:(M),each=n.all/2),rep(M+1,times=n.control/2),rep(M+2,times=n.recovered/2))
      A[,"Piece"] = c(rep(1:(n.all/2),times=M),rep(1:(n.control/2),times=1),rep(1:(n.recovered/2),times=1))
      
      meanCentred_test = cbind(A[,c(1:3)],A[,-c(1:3)]-matrix(center,nrow=nrow(A),ncol=ncol(A)-3,byrow=T))
      projections = as.matrix(meanCentred_test[,-c(1:3)])%*%eigen_vectors ##etap 6
      P = cbind(A[,1:3],projections) ##etap 6
      
      vars = 1
      
      LR.final.Nor=LR.final.KDE = 1
      
      for (v in 1:length(vars))
      {
        variables = vars[v]+3
        p = length(variables) 
        
        y.1 = P[which(P$Item==M+2),]
        y.2 = P[which(P$Item==M+1),]
        
        population = P[which(!P$Item %in% c(M+1,M+2)),]  
        
        n = length(unique(population$Piece))
        m = length(unique(population$Item))
        
        results.BW = BW(population=population,variables=variables,n=n,m=m) 
        W = results.BW$W
        B = results.BW$B
        
        ifelse(p>1,x.mean<-colMeans(population[,variables]),x.mean<-mean(population[,variables]))
        
        h = (4/(m*(2*p+1)))^(1/(p+4))
        
        k.1 = length(y.1$Item) 
        k.2 = length(y.2$Item) 
        
        y.mean.1 = matrix(apply(as.matrix(y.1[,variables]), 2, mean), nrow = 1) 
        y.mean.2 = matrix(apply(as.matrix(y.2[,variables]), 2, mean), nrow = 1) 
        y.star = (k.1*y.mean.1+k.2*y.mean.2)/(k.1+k.2)    
        
        results.LR.Nor = LR.Nor.function(y.mean.1, y.mean.2, y.star, x.mean, W, B, k.1, k.2) ##etap 7
        LR.Nor = results.LR.Nor$LR.Nor
        
        results.LR.KDE= LR.KDE.function(y.mean.1, y.mean.2, y.star, W, B, h, population, variables, p, k.1, k.2) ##etap 7
        LR.KDE = results.LR.KDE$LR.KDE
        
        LR.final.Nor=LR.final.Nor*LR.Nor
        LR.final.KDE=LR.final.KDE*LR.KDE
      }
      output.matrix.Nor[j,i,s] = LR.final.Nor
      output.matrix.KDE[j,i,s] = LR.final.KDE
    }
  }
}

saveRDS(output.matrix.Nor,file="hybridLR_comparison_casework_Nor.txt")
saveRDS(output.matrix.KDE,file="hybridLR_comparison_casework_KDE.txt")

