#Author: Penelope A. Hancock
library(INLA)

inla_marker_frequency<-function(data1){
	# data1 is an R data frame containing the data in one of the files kdr_frequenciesF_gamb.csv, 
	# Cyp6p4_236M_frequencies_gamb.csv or Cyp4j5_43F_frequencies_gamb.csv (all can be found in 
	# the folder "allele_frequencies".
	# The version of these files in this repository have had the latitudes and longitudes obfuscated
	# for data protection.

	long_lat_all<-cbind(data1[,"lon"],data1[,"lat"])

	#INLA mesh parameters
	m1.cutoff<-0.2
	m1.min.angle<-c(25,25)
	m1.max.edge<-c(0.05,500)
	tmesh.st<-1
	tmesh.end<-5
	tmesh.end2<-5
	tmesh.by<-1

	#Make the spatial mesh
	mesh = inla.mesh.2d(loc=cbind(long_lat_all[,1],long_lat_all[,2]),
						cutoff=m1.cutoff,
						min.angle=m1.min.angle,
						max.edge=m1.max.edge)

	#Make the time mesh
	mesh1d=inla.mesh.1d(seq(tmesh.st,tmesh.end,by=tmesh.by),interval=c(tmesh.st,tmesh.end2),degree=2, boundary='free')

	#Create the SPDE model using penalized complexity priors
	IRT_spde<- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
								   prior.range = c(0.1,0.1), 
								   prior.sigma = c(5,0.1))

	IRT_dat1 <- data.frame(y=as.vector(round(data1[,"total_freq"]*2*data1[,"tot_mosq"])),n=2*data1[,"tot_mosq"], w=rep(1,length(data1[,1])), rnd = data1[,"rnd"], xcoo=long_lat_all[,1],
						   ycoo=long_lat_all[,2])

	#Make an index for the spatial Gaussian Markov random effect parameters
	IRT_iset1 <- inla.spde.make.index('one.field', n.spde=mesh$n, n.group=mesh1d$m)

	#Make the projector matrix
	IRT_A1 <- inla.spde.make.A(mesh,loc=cbind(IRT_dat1$xcoo, IRT_dat1$ycoo),group=IRT_dat1$rnd,group.mesh=mesh1d)

	#Combine everything into an INLA stack object
	IRT_stk1 <- inla.stack(tag='one.data', data=list(y=IRT_dat1$y,n=IRT_dat1$n), A=list(IRT_A1,1),effects=list(IRT_iset1, b0.1=IRT_dat1$w)) 

	#Set prior parameters for the AR1 temporally autocorrelated random effect
	ar1_prior_pars<- list(theta=list(initial=0.3, param=c(0, 1)))

	bin_form <- y ~ 0 + b0.1 + f(one.field, model=IRT_spde, group=one.field.group,control.group=list(model='ar1',hyper=ar1_prior_pars)) 
	bin.res <- inla(bin_form, family=c('Binomial'),Ntrials = n, data=inla.stack.data(IRT_stk1), control.predictor=list(compute=TRUE,A=inla.stack.A(IRT_stk1)),control.compute=list(config=TRUE),control.inla=list(int.strategy='eb'))

	return(bin.res)

}

