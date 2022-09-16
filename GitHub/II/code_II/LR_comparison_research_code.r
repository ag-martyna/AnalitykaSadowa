source("../BW_function.R")
source("../LR_Nor_comparison_function.R")
source("../LR_KDE_comparison_function.R")

if (!require("mvtnorm")) install.packages("mvtnorm")

population = read.table("../demo_glass_database.txt", header = TRUE)
data.analysed = population

items = unique(population$Item)
M = length(unique(population$Item))
n = length(unique(population$Piece))

variables.list = list(c("logNaO","logSiO","logCaO"),c("logMgO")) 
ND = c("N","D")
variable.name = "NaSiCaMg"

output.matrix.Nor = output.matrix.KDE = matrix(1,ncol = M, nrow = M,dimnames=list(c(as.character(unique(population$Name))),c( as.character(unique(population$Name)))))

for (v in 1:length(variables.list))
{
  variables = variables.list[[v]]
  variables = which(colnames(population) %in% variables)
  p = length(variables) 
  
  for (i in 1:M)
  {  
    for (j in i:M) 
    {	
      if (i == j) 
      {
        y.1.2 = data.analysed[which(data.analysed$Item == items[i]),] 
        y.1 = data.frame(y.1.2[c(1,3,5,7,9,11),]) 
        y.2 = data.frame(y.1.2[c(2,4,6,8,10,12),]) 
        
        population = data.analysed[which(data.analysed$Item != items[i]),] ##A, gdy H1 jest prawdziwa
        m = length(unique(population$Item)) 
      }
      else 
      {
        y.1 = data.frame(data.analysed[which(data.analysed$Item == items[i]),]) 
        y.2 = data.frame(data.analysed[which(data.analysed$Item == items[j]),]) 
        
        population = data.analysed[which(!data.analysed$Item %in% c(items[i],items[j])),] ##A, gdy H2 jest prawdziwa
        m = length(unique(population$Item)) 
      }
      
      results.BW = BW(population, variables, n, m) 
      W = results.BW$W
      B = results.BW$B
            
      ifelse(p>1,x.mean<-colMeans(population[,variables]),x.mean<-mean(population[,variables]))
      
      h = (4/(m*(2*p+1)))^(1/(p+4))
       
      k.1 = length(y.1$Item) 
      k.2 = length(y.2$Item) 
      
      y.mean.1 = matrix(apply(as.matrix(y.1[,variables]), 2, mean), nrow = 1) 
      y.mean.2 = matrix(apply(as.matrix(y.2[,variables]), 2, mean), nrow = 1) 
      y.star = (k.1*y.mean.1+k.2*y.mean.2)/(k.1+k.2)    
      
      results.LR.Nor = LR.Nor.function(y.mean.1, y.mean.2, y.star, x.mean, W, B, k.1, k.2)
      LR.Nor = results.LR.Nor$LR.Nor
      
      results.LR.KDE= LR.KDE.function(y.mean.1, y.mean.2, y.star, W, B, h, population, variables, p, k.1, k.2)
      LR.KDE = results.LR.KDE$LR.KDE
      
      if (ND[v] == "N") 
      {
        output.matrix.Nor[i,j] = output.matrix.Nor[i,j] * LR.Nor
        output.matrix.KDE[i,j] = output.matrix.KDE[i,j] * LR.KDE
      } else
      {
        output.matrix.Nor[i,j] = output.matrix.Nor[i,j] * 1/LR.Nor
        output.matrix.KDE[i,j] = output.matrix.KDE[i,j] * 1/LR.KDE
      } 
    }	
  }	
}

write.table(signif(output.matrix.Nor, digits = 4), file = paste(variable.name,"_comparison_research_Nor.txt", sep=""), quote = FALSE, sep = "\t")
write.table(signif(output.matrix.KDE, digits = 4), file = paste(variable.name,"_comparison_research_KDE.txt", sep=""), quote = FALSE, sep = "\t")

fp.all = M*(M-1)/2 
LR.H2.Nor = output.matrix.Nor[upper.tri(output.matrix.Nor, diag = FALSE)]
LR.H2.KDE = output.matrix.KDE[upper.tri(output.matrix.KDE, diag = FALSE)]
LR.H1.Nor = diag(output.matrix.Nor) 
LR.H1.KDE = diag(output.matrix.KDE) 
fp.Nor = length(which(LR.H2.Nor>1))/fp.all*100
fp.KDE = length(which(LR.H2.KDE>1))/fp.all*100
fn.Nor = length(which(LR.H1.Nor<1))/M*100
fn.KDE = length(which(LR.H1.KDE<1))/M*100

error_rates = matrix(0, nrow = 4, ncol = 1,dimnames=list(c("fp_Nor", "fp_KDE", "fn_Nor", "fn_KDE"),c(variable.name)))
error_rates[1,1] = fp.Nor
error_rates[2,1] = fp.KDE
error_rates[3,1] = fn.Nor
error_rates[4,1] = fn.KDE
write.table(signif(error_rates, digits = 3), file=paste(variable.name,"_comparison_research_error_rate.txt",sep=""), quote = FALSE, sep = "\t")

source("../ECE_function.R")
ECE = ECE_plot(LR.H1.Nor,LR.H2.Nor)
dev.copy(png, paste(variable.name,"_Nor_ECE.png", sep=""))
dev.off()

ECE = ECE_plot(LR.H1.KDE,LR.H2.KDE)
dev.copy(png, paste(variable.name,"_KDE_ECE.png", sep=""))
dev.off()



