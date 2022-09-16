BW = function(population,variables,n,m)
{
	if (n==1) W=0 else
	W = Reduce("+",lapply(split(population[,variables], population$Item), function(X) var(X)))

	W = W/m

	means = aggregate(. ~ population$Item, as.data.frame(population[,variables]), mean)
	B = var(means[,-1])-W/n

	return (list(W = W, B = B))
}