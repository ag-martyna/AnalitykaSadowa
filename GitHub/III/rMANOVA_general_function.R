rMANOVA_general_function = function(data)
{
	if (!require("rARPACK")) install.packages("rARPACK")

	rMANOVA_function = function(database,variables)
	{
		database = database[,c(1:3,variables)]
		m=length(unique(database$Item))
		N=nrow(database)
		p=length(variables)
		n =N/m
		sample_means = aggregate(. ~ Item, database[,-c(1,3)] , mean)[,-1]
		B=n*var(sample_means)*(m-1) 
		W=matrix(0,nrow=p,ncol=p)
		for (ii in 1:length(unique(database$Item)))
		{
			object = database[which(database$Item==unique(database$Item)[ii]),]
			W = W+var(object[,-c(1:3)])*(n-1)
		}

		sample_means_matrix = apply(sample_means,2,function(x) rep(x,each=n))
		A_centered = database[,-c(1:3)]-sample_means_matrix

		var_wij = matrix(NA,p,p)
		for (iii in 1:(p))
		{
			if (iii==p) var_wij[iii,iii:p]  = var(A_centered[,iii:p]*A_centered[,iii:p]) else
			var_wij[iii,iii:p]  = apply(apply(A_centered[,iii:p],2,function(x) A_centered[,iii]*x),2,var)
		}
		var_wij=N*var_wij 

		T = diag(diag(W),nrow=p,ncol=p)
		delta = sum(var_wij[col(var_wij)>row(var_wij)])/sum(W[col(W)>row(W)]^2)
		delta = max(0,min(delta,1))

		Jstar = tryCatch(solve((1-delta)*W+delta*T)%*%B, error = function(e) {
		delta=1;print(delta)
		solve((1-delta)*W+delta*T)%*%B})

		eigen = eigs(Jstar,k=10)
		eigen_vectors = eigen$vectors
		eigen_values = eigen$values

		return(list(eigen_vectors=eigen_vectors,eigen_values=eigen_values))
	}

	input.data = data
	train=input.data

	center = apply(train[,-c(1:3)],2,mean)
	meanCentred = cbind(train[,c(1:3)],train[,-c(1:3)]-matrix(center,nrow=nrow(train),ncol=ncol(train)-3,byrow=T))
	rMANOVA_res = tryCatch(rMANOVA_function(database=meanCentred,variables=4:ncol(input.data)), error = function(e) 
	{b0=length(which(input.data[,-c(1:3)]==0))
	input.data[,-c(1:3)][(input.data[,-c(1:3)]==0)]=rnorm(b0,mean=10^-8,sd=10^-8)
	train=input.data
	center = apply(train[,-c(1:3)],2,mean)
	meanCentred = cbind(train[,c(1:3)],train[,-c(1:3)]-matrix(center,nrow=nrow(train),ncol=ncol(train)-3,byrow=T))
	rMANOVA_function(database=meanCentred,variables=4:ncol(input.data))})

	eigen_vectors = Re(rMANOVA_res$eigen_vectors)
	eigen_values = Re(rMANOVA_res$eigen_values)

	list(results = eigen_vectors,center = center)
}



