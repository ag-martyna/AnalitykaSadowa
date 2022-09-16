LR.KDE.function = function(y.mean.1, y.mean.2, y.star, W, B, h, population, variables, p, k.1, k.2)
{
	density.1.1 = mvtnorm:::dmvnorm(y.mean.1,mean=y.mean.2,sigma=as.matrix(W/k.1+W/k.2))

	means = aggregate(.~population$Item,as.data.frame(population[,variables]),mean)[,-1] 
	if (p==1)  means = matrix(means,ncol=1)

	density.1.2 = mean(apply(means,1,function(x) mvtnorm:::dmvnorm(y.star,mean=x,sigma=as.matrix(W/(k.1+k.2)+B*h^2))))  

	density.2.1 = mean(apply(means,1,function(x) mvtnorm:::dmvnorm(y.mean.1,mean=x,sigma=as.matrix(W/k.1+B*h^2))))  
	density.2.2 = mean(apply(means,1,function(x) mvtnorm:::dmvnorm(y.mean.2,mean=x,sigma=as.matrix(W/k.2+B*h^2))))  

	LR.KDE = density.1.1*density.1.2/density.2.1/density.2.2

	result = list(LR.KDE = LR.KDE)
	return (result)
}























