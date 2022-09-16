source("../BW_function.R")
source("../LR_Nor_comparison_function.R")
source("../LR_KDE_comparison_function.R")

if (!require("mvtnorm")) install.packages("mvtnorm")


population = read.table("../demo_glass_database.txt", header = TRUE) ##A
data.recovered = read.table("../data_recovered.txt", header = TRUE) 
data.control = read.table("../data_control.txt", header = TRUE) 

M = length(unique(population$Item)) 
m.control = length(unique(data.control$Item)) 
m.recovered = length(unique(data.recovered$Item)) 

n = length(unique(population$Piece))

variables.list = list(c("logNaO","logSiO","logCaO"),c("logMgO")) 
ND = c("N","D")
variable.name = "NaSiCaMg"

output.matrix.Nor = output.matrix.KDE = matrix(1, ncol = m.control, nrow = m.recovered,dimnames=list(c(as.character(unique(data.recovered$Name))),c(as.character(unique(data.control$Name)))))

for (v in 1:length(variables.list))
{
  variables = variables.list[[v]]
  variables = which(colnames(population) %in% variables)
  p = length(variables)
  
  results.BW = BW(population, variables, n, M) 
  W = results.BW$W
  B = results.BW$B

  ifelse(p>1,x.mean<-colMeans(population[,variables]),x.mean<-mean(population[,variables]))
  
  h = (4/(m*(2*p+1)))^(1/(p+4))
  
  for (i in 1:m.recovered)
  {
  	y.1 = data.frame(data.recovered[which(data.recovered$Item == i),]) 
  
  	for (j in 1:m.control)
  	{		
  		y.2 = data.frame(data.control[which(data.control$Item == j),]) 
  
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

write.table(signif(output.matrix.Nor, digits = 4), file = paste(variable.name,"_comparison_casework_Nor.txt", sep=""), quote = FALSE, sep = "\t")
write.table(signif(output.matrix.KDE, digits = 4), file = paste(variable.name,"_comparison_caseworkKDE.txt", sep=""), quote = FALSE, sep = "\t")



