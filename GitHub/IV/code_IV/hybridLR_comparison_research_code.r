source("../BW_function.R")
source("../LR_Nor_comparison_function.R")
source("../LR_KDE_comparison_function.R")
source("../rMANOVA_general_function.R")

if (!require("mvtnorm")) install.packages("mvtnorm")


population = read.table(paste("../spectra_database.txt",sep=""),header=T)
data.analysed = population

M = length(unique(population$Name))
n.all=length(unique(population$Piece))

validationSetsNo = 5 


output.matrix.Nor = output.matrix.KDE = array(NA,c(M,M,validationSetsNo),dimnames=list(as.character(unique(population$Name)),as.character(unique(population$Name)),1:validationSetsNo))

for (s in 1:validationSetsNo)
{
  tvec = NULL
  n.val = n.all/3
  for (g in 1:M)
  {tvec = c(tvec,sample(rep(c(1,2,3),each=n.val)))}
  
  A.1 = data.analysed[which(tvec == 1),] ##etap 1 na Rys. 1.3 i 1.4
  A.1[,"Piece"] = rep(1:n.val,times=M)
  A.2 = data.analysed[which(tvec == 2),] ##etap 1 na Rys. 1.3 i 1.4
  A.2[,"Piece"] = rep(1:n.val,times=M)
  A.3 = data.analysed[which(tvec == 3),] ##etap 1 na Rys. 1.3 i 1.4
  A.3[,"Piece"] = rep(1:n.val,times=M)
  
  rMANOVA_res = rMANOVA_general_function(data=as.data.frame(A.1)) ##etap 2 na Rys. 1.4
  eigen_vectors_ij = rMANOVA_res$results
  center_ij = rMANOVA_res$center
  
  for (i in 1:M)
  {  
    
    ##jeœli j/=i
    B = A.1[-which(A.1$Item == i),] ##etap 2 na Rys. 1.3
    B$Item = rep(1:(M-1),each=n.val)
    
    rMANOVA_res = rMANOVA_general_function(data=as.data.frame(B)) ##etap 3 na Rys. 1.3
    eigen_vectors_inotj = rMANOVA_res$results
    center_inotj = rMANOVA_res$center
    
    meanCentred_test_inotj = cbind(A.2[,c(1:3)],A.2[,-c(1:3)]-matrix(center_inotj,nrow=nrow(A.2),ncol=ncol(A.2)-3,byrow=T))
    projections_inotj = as.matrix(meanCentred_test_inotj[,-c(1:3)])%*%eigen_vectors_inotj ##etap 4 na Rys. 1.3
    
    P_inotj = cbind(A.2[,1:3],projections_inotj) ##etap 4 na Rys. 1.3
    
    ##jeœli j==i
    A = rbind(A.2,A.3[which(A.3$Item==i),]) ##etap 3 na Rys. 1.4
    A[,"Name"] = c(A.2[,3],rep(1000,times=n.val))
    A[,"Item"] = rep(1:(M+1),each=n.val)
    A[,"Piece"] = rep(1:n.val,times=M+1)
    
    meanCentred_test_ij = cbind(A[,c(1:3)],A[,-c(1:3)]-matrix(center_ij,nrow=nrow(A),ncol=ncol(A)-3,byrow=T))
    projections_ij = as.matrix(meanCentred_test_ij[,-c(1:3)])%*%eigen_vectors_ij ##etap 4 na Rys. 1.4
    
    P_ij = cbind(A[,1:3],projections_ij) ##etap 4 na Rys. 1.4		
    
    for (j in 1:M) 
    {  
      if (i==j) P = P_ij else P = P_inotj
      
      vars = 1
      
      LR.final.Nor=LR.final.KDE = 1
      
      for (v in 1:length(vars))
      {
        variables = vars[v]+3
        p = length(variables) 
        
        if (i==j)
        {
          y.1 = P[which(P$Item == (M+1)),]
          y.2 = P[which(P$Item == j),]
          
          population = P[which(!P$Item %in% c(j,M+1)),]    
        } else
        {
          y.1 = P[which(P$Item==i),]
          y.2 = P[which(P$Item==j),]
          
          population = P[which(!P$Item %in% c(i,j)),]  
        }
        
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
        
        results.LR.Nor = LR.Nor.function(y.mean.1, y.mean.2, y.star, x.mean, W, B, k.1, k.2) ##etap 5 na Rys. 1.3 i 1.4
        LR.Nor = results.LR.Nor$LR.Nor
        
        results.LR.KDE= LR.KDE.function(y.mean.1, y.mean.2, y.star, W, B, h, population, variables, p, k.1, k.2) ##etap 5 na Rys. 1.3 i 1.4
        LR.KDE = results.LR.KDE$LR.KDE
        
        LR.final.Nor=LR.final.Nor*LR.Nor
        LR.final.KDE=LR.final.KDE*LR.KDE
      }
      
      output.matrix.Nor[i,j,s] = LR.final.Nor
      output.matrix.KDE[i,j,s] = LR.final.KDE
    }
  }
}

saveRDS(output.matrix.Nor,file="hybridLR_comparison_research_Nor.txt")
saveRDS(output.matrix.KDE,file="hybridLR_comparison_research_KDE.txt")






error_rates = matrix(0, nrow = 4, ncol = validationSetsNo,dimnames=list(c("fp_Nor", "fp_KDE", "fn_Nor", "fn_KDE"),1:validationSetsNo))
source("../ECE_function.R")

for (ss in 1:validationSetsNo)
{
  true.H1.Nor = diag(output.matrix.Nor[,,ss])
  true.H2.Nor = output.matrix.Nor[,,ss][row(output.matrix.Nor[,,ss])!=col(output.matrix.Nor[,,ss])]
  true.H1.KDE = diag(output.matrix.KDE[,,ss])
  true.H2.KDE = output.matrix.KDE[,,ss][row(output.matrix.KDE[,,ss])!=col(output.matrix.KDE[,,ss])]
  fn.Nor = length(which(true.H1.Nor<1))/M*100
  fp.Nor = length(which(true.H2.Nor>1))/(M*(M-1))*100
  fn.KDE = length(which(true.H1.KDE<1))/M*100
  fp.KDE = length(which(true.H2.KDE>1))/(M*(M-1))*100
  
  error_rates[1,ss] = fp.Nor
  error_rates[2,ss] = fp.KDE
  error_rates[3,ss] = fn.Nor
  error_rates[4,ss] = fn.KDE
  
  write.table(signif(error_rates, digits = 3), file="hybridLR_comparison_research_error_rate.txt", quote = FALSE, sep = "\t")
  
  ECE_results = ECE_plot(true.H1.Nor,true.H2.Nor)
  dev.copy(png,paste("hybridLR_comparison_research_Nor_ECE_s",ss,".png",sep=""));dev.off()
  
  ECE_results = ECE_plot(true.H1.KDE,true.H2.KDE)
  dev.copy(png,paste("hybridLR_comparison_research_KDE_ECE_s",ss,".png",sep=""));dev.off()
}

