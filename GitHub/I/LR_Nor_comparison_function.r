LR.Nor.function = function(y.mean.1, y.mean.2, y.star, x.mean, W, B, k.1, k.2)
{
	density.1.1 = mvtnorm:::dmvnorm(y.mean.1,mean=y.mean.2,sigma=as.matrix(W/k.1+W/k.2))
	density.1.2 = mvtnorm:::dmvnorm(y.star,mean=x.mean,sigma=as.matrix(W/(k.1+k.2)+B))
	
	density.2.1 = mvtnorm:::dmvnorm(y.mean.1,mean=x.mean,sigma=as.matrix(W/k.1+B))
	density.2.2 = mvtnorm:::dmvnorm(y.mean.2,mean=x.mean,sigma=as.matrix(W/k.2+B))
	 
	LR.Nor = density.1.1*density.1.2/density.2.1/density.2.2

	result = list(LR.Nor = LR.Nor)
	return (result)
}

