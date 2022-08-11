mod <- function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		
		
		if(t - timeg < 1.5 ){
			CSFY = 1
		}else{
			CSFY = 0
		}
		
		
		# #DELTIME = TIME-TDOS
		# IF(t < 1.5) {
		#   CSFY = 1
		# }
		# 
		
		ID50C = ID50P
		
		#print(state)
		
		DA_1DT = -KDEP*A_1
		DA_2DT = -KDEC*A_2
		VIRP    = KDEP*A_1
		VIRC    = KDEC*A_2
		PCEFF = 1-(VIRP/ID50P+VIRC/ID50C)/(1+VIRP/ID50P+VIRC/ID50C) #Assumes U50=1 in Minto model
		
		
		
		
		KTR  = 4/MTT*(1+INCTRANS*CSFY)
		KCIRC = 4/MTT
		DA_3DT = KTR*A_3*(1+INCPROL*CSFY)*PCEFF*(BASE/A_7)**POWER-KTR*A_3
		DA_4DT = KTR*A_3-KTR*A_4
		DA_5DT = KTR*A_4-KTR*A_5
		DA_6DT = KTR*A_5-KTR*A_6
		DA_7DT = KTR*A_6- KCIRC*A_7
		
		
		list(c(DA_1DT, DA_2DT, DA_3DT, DA_4DT, DA_5DT, DA_6DT, DA_7DT))
		
	})
}


sim <- function(state, totaldose,parameters,maxtime) {
	

	
	out <- data.frame()
	# parameters <- c(BASE = BASE,
	#                 MTT = MTT,
	#                 POWER = POWER,
	#                 ID50P = ID50P,
	#                 KDEP = KDEP,
	#                 KDEC = KDEC,
	#                 INCTRANS = INCTRANS,
	#                 INCPROL =INCPROL)

	
	# state <- c(A_1 = DOSEP, 
	#            A_2 = 0, 
	#            A_3 = BASE,
	#            A_4 = BASE,
	#            A_5 = BASE,
	#            A_6 = BASE,
	#            A_7 = BASE)
	
	time.lt.maxtime <- totaldose$time[totaldose$time < maxtime]
	
	if(totaldose$time[1] != 0){
		out <- rbind(out,
								 data.frame(seq(0, totaldose$time[1], by=0.5),
								 					 parameters[1],
								 					 parameters[1],
								 					 parameters[1],
								 					 parameters[1],
								 					 parameters[1],
								 					 parameters[1],
								 					 parameters[1]
								 ))
		colnames(out) <-  c("time", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7")
	}
	

	
	timeg = -100000
	for (i in 1:length(time.lt.maxtime)) {
		
			
		
		if(totaldose$time[i] != 0){
			state <- as.numeric(out[nrow(out),-1])
			names(state) <- c("A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7")
		}
		
		if(totaldose$variable[i] == "dosep"){
			state["A_1"] <- state["A_1"] + totaldose$dose[i]
		}else if (totaldose$variable[i] == "dosec"){
			state["A_2"] <- state["A_2"] + totaldose$dose[i]
		}else if (totaldose$variable[i] == "doseg"){
			timeg = totaldose$time[i]
		}
		
		if(i !=length(time.lt.maxtime)){
			times <- seq(totaldose$time[i], totaldose$time[i+1], by=0.5)
		}else{
			times <- seq(totaldose$time[i], maxtime, by=0.5)
		}
		
		#print(state)    
		
		
		parameters.sim <- c(parameters, timeg= timeg) 
		#print(state)
		df <- as.data.frame(ode(y = state, times = times, func = mod, parms = parameters.sim))
		out <- rbind(out, df)
	} 
	return(out)
	
}

# Newton algorithm
newton <- function(f, tol = 1e-7, x0, N=30) {#1e-7
	h = 1e-7
	i = 1; x1 = x0
	p = matrix(0, N +10, length(x0))
	colnames(p) <- names(x0)
	while (i <= N + 10) {
		G = jacobian(f, x0)
		H = hessian(f, x0)
		H_inv <- solve(H)
		x1 = x0 - H_inv%*%t(G)
		p[i,] = t(x1)
		i = i + 1
		print(i)
		#print(norm(as.matrix(x1 - x0)))
		#print(solve(H)%*%t(G))
		print(G)
		print(exp(x1))
		if (norm(as.matrix(x1 - x0)) < tol) break
		
		
		if(i > N){
			H <- hessian(f, x1)
			G <- jacobian(f, x1)
			H_inv <- solve(H)
			if(is.positive.definite(1/2*(H_inv+t(H_inv)))) break
		} 
		x0 = x1
	}
	
	return(list(p = p[1:(i-1),],
							G = jacobian(f, x1),
							H = hessian(f, x1)))
}

runmain <- function(totaldose,  SEX, DM,  maxtime, BASE, DV = NULL, plot.variability ) {
	
	
	# Number of individuals: 173
	# 
	# THETA                               OMEGA              SIGMA    
	# T1. BASE 5.124     5.324  (0.02832)         O1. BASE  0.2374   (0.115)             
	# T2. MTT 5.4D      4.668    (0.031)          O2. MTT  0.1608  (0.1533)             
	# T3. POWER 0.23    0.1953  (0.06768)  O3. ID50P=ID50C  0.2722  (0.1899)             
	# T4. ID50P 7.17     88.48  (0.08119)         O4. KDEP  0.8246  (0.1565)             
	# T5. KDEP          0.03704  (0.08196)         O5. KDEC  0.9515  (0.2367)             
	# T6. KDEC           0.1718   (0.2328)     O6. INCTRANS  0.3388   (0.491)             
	# T7. Increase of KTR 3.073   (0.1716)                                                
	# T8. Increase of Kprol 0.208   (0.1569)                                                
	# T9. PROP   0.2941   (0.04165)                                                
	# T10. ADD.    0.529  (0.08584)                                                
	# ID50PSEX1           -0.3338   (0.1766)                                                
	# ID50PDM1             0.4596   (0.3408)                                               
	# 
	
	#without GCSF IIV
	# THETA                               OMEGA              SIGMA    
	# T1. BASE 5.124    5.341  (0.02815)         O1. BASE  0.2301  (0.1142)             
	# T2. MTT 5.4D    4.639  (0.03299)          O2. MTT  0.1616  (0.1632)             
	# T3. POWER 0.23   0.1882  (0.06447)  O3. ID50P=ID50C  0.2332  (0.2454)             
	# T4. ID50P 7.17    88.88  (0.07822)         O4. KDEP  0.9557  (0.1378)             
	# T5. KDEP  0.03261    (0.132)         O7. KDEC    0.92  (0.2409)             
	# T6. KDEC   0.1877   (0.2451)                                                
	# T7. Increase of KTR    3.502  (0.08684)                                                
	# T8. Increase of Kprol   0.2172  (0.08525)                                                
	# T9. PROP   0.2998  (0.03937)                                                
	# T10. ADD.   0.5271  (0.08519)                                                
	# ID50PSEX1  -0.3338    (0.164)                                                
	# ID50PDM1   0.4849    (0.315)                       
	
	GCSFIIV = 0
	
	if(GCSFIIV){
		TVBASE <- 5.324
		BASE <- BASE
		TVMTT <- 4.668
		TVPOWER  <-  0.1953
		TVID50P <-  88.48
		TVKDEP <- 0.03704
		TVKDEC <- 0.1718
		TVINCTRANS <- 3.073
		TVINCPROL <- 0.208
		ID50P <- TVID50P*(1+0.4596*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2374,  MTT=0.1608, ID50P =0.2722, KDEP=0.8246,KDEC= 0.9515,INCTRANS= 0.3388,INCPROL= 0.3388)
		sigadd <- 0.529
		sigprop <- 0.2941
		
	}else{
		TVBASE <- 5.341
		BASE <- BASE
		TVMTT <- 4.639
		TVPOWER  <-  0.1882
		TVID50P <-  88.88
		TVKDEP <- 0.03261
		TVKDEC <- 0.1877
		TVINCTRANS <- 3.502
		TVINCPROL <- 0.2172
		ID50P <- TVID50P*(1+0.4849*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2301, MTT= 0.1616,ID50P = 0.2332,KDEP= 0.9557,  KDEC= 0.92)
		sigadd <- 0.5271
		sigprop <- 0.2998
		
	}
	
	
	######################################
	#####population mean simulation ######
	######################################
	
	parameters <- c(BASE = BASE,
									MTT = TVMTT,
									POWER = TVPOWER,
									ID50P = ID50P,
									KDEP = TVKDEP,
									KDEC = TVKDEC,
									INCTRANS = TVINCTRANS,
									INCPROL =TVINCPROL)
	
	state <- c(A_1 = 0, 
						 A_2 = 0, 
						 A_3 = BASE,
						 A_4 = BASE,
						 A_5 = BASE,
						 A_6 = BASE,
						 A_7 = BASE)
	
	
	out <- sim(state,totaldose, parameters, maxtime)
	#############################
	###### variabilbity #########
	#############################
	
	if(as.numeric(plot.variability) == 1){
		
		Nrand <- 500
		set.seed(1111)
		
		if(GCSFIIV){
			prm.eta <- data.frame(ETA_BASE = rep(NA,Nrand),
														ETA_MTT = rep(NA,Nrand),
														ETA_ID50 = rep(NA,Nrand),
														ETA_KDEP = rep(NA,Nrand),
														ETA_KDEC = rep(NA,Nrand),
														ETA_TRANS = rep(NA,Nrand),
														ETA_PROL = rep(NA,Nrand))
			for(i in 1:7){  
				
				prm.eta[,i ] <- parameters[-3][i]*exp(rnorm(Nrand,0, omega[[i]]))
			}
			
		}else{
			prm.eta <- data.frame(ETA_BASE = rep(NA,Nrand),
														ETA_MTT = rep(NA,Nrand),
														ETA_ID50 = rep(NA,Nrand),
														ETA_KDEP = rep(NA,Nrand),
														ETA_KDEC = rep(NA,Nrand))
			for(i in 1:5){  
				prm.eta[,i ] <- parameters[-3][i]*exp(rnorm(Nrand,0, omega[[i]]))
			}
		}
	
		
		###simulation for interindividual variability
		for(idx.var in 1:Nrand){
			if(GCSFIIV){
				parameters.eta <- c(BASE = BASE,
														MTT = prm.eta[idx.var,2],
														POWER = TVPOWER,
														ID50P = prm.eta[idx.var,3],
														KDEP = prm.eta[idx.var,4],
														KDEC = prm.eta[idx.var,5],
														INCTRANS = prm.eta[idx.var,6],
														INCPROL =prm.eta[idx.var,7])
			}else{
				parameters.eta <- c(BASE = BASE,
														MTT = prm.eta[idx.var,2],
														POWER = TVPOWER,
														ID50P = prm.eta[idx.var,3],
														KDEP = prm.eta[idx.var,4],
														KDEC = prm.eta[idx.var,5],
														INCTRANS = TVINCTRANS,
														INCPROL =TVINCPROL)
			}
		
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = BASE,
								 A_4 = BASE,
								 A_5 = BASE,
								 A_6 = BASE,
								 A_7 = BASE)
			
			
		
			temp.out <- sim(state,totaldose, parameters.eta, maxtime)
			out <- cbind(out, var = temp.out[,"A_7"])
		}
		colnames(out)[9:508] <- paste0("var", 1:500)
	}
	
	
	
	###############################################
	#######individualized parameter estimation#####
	###############################################
	
	if(!is.null(DV)){

		parameters <- c(BASE = TVBASE,
										MTT = TVMTT,
										POWER = TVPOWER,
										ID50P = ID50P,
										KDEP = TVKDEP,
										KDEC = TVKDEC,
										INCTRANS = TVINCTRANS,
										INCPROL =TVINCPROL)
		
		ETA_BASE <- 0
		ETA_MTT <- 0
		ETA_ID50P <- 0
		ETA_KDEP <- 0
		ETA_KDEC <- 0
		ETA_INCTRANS <- 0
		ETA_INCPROL <- 0
		
		# Negative Log-likelihood
		llmod <- function(state,totaldose,parameters, maxtime){
			ind_pred <- sim(state,totaldose,parameters, maxtime)
			ind_pred <- ind_pred %>% filter(time %in% DV$time) %>% select(time, A_7) %>% mutate(y = A_7) %>% select(-A_7) %>% unique()
			like <- -dnorm((DV$y-ind_pred$y)/sqrt(sigprop^2*(ind_pred$y)^2+sigadd^2), 0, 1, log=T)
			#SQRT(THETA(9)**2*IPRED**2+THETA(10)**2)
			
			return(like)
		}
		
		ll <- function(parameters) {
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = parameters[["BASE"]],
								 A_4 = parameters[["BASE"]],
								 A_5 = parameters[["BASE"]],
								 A_6 = parameters[["BASE"]],
								 A_7 = parameters[["BASE"]])
			temp.sum <- sum(llmod(state,totaldose,parameters, max(DV$time))) - 
				dnorm(log(parameters["BASE"]), log(TVBASE), omega["BASE"], log=T) -
				dnorm(log(parameters["MTT"]), log(TVMTT), omega["MTT"], log=T) -
				dnorm(log(parameters["ID50P"]), log(ID50P), omega["ID50P"], log=T) -
				dnorm(log(parameters["KDEP"]), log(TVKDEP), omega["KDEP"], log=T) -
				dnorm(log(parameters["KDEC"]), log(TVKDEC), omega["KDEC"], log=T) -
				dnorm(log(parameters["INCTRANS"]), log(TVINCTRANS), omega["INCTRANS"], log=T) -
				dnorm(log(parameters["INCPROL"]), log(TVINCPROL), omega["INCPROL"], log=T)
			return(temp.sum)
		}
		
		ll2 <- function(x){
			if(GCSFIIV){
				parameters <- c(BASE = exp(x[1]),
												MTT = exp(x[2]),
												POWER = TVPOWER,
												ID50P = exp(x[3]),
												KDEP = exp(x[4]),
												KDEC = exp(x[5]),
												INCTRANS = exp(x[6]),
												INCPROL = exp(x[7]))
			}else{
				parameters <- c(BASE = exp(x[1]),
												MTT = exp(x[2]),
												POWER = TVPOWER,
												ID50P = exp(x[3]),
												KDEP = exp(x[4]),
												KDEC = exp(x[5]),
												INCTRANS = TVINCTRANS,
												INCPROL = TVINCPROL)
			}
			
			return(ll(parameters))
		}
		#ll2(log(c(5.324, 4.668, 88.48, 0.03704, 0.1718, 3.073, 0.208)))
		
		
		#newton.results <- newton(ll2, x0=log(c(5.324, 4.668, ID50P, 0.03704, 0.1718, 3.073, 0.208)))
		if(GCSFIIV){
			newton.results <- newton(ll2, x0=log(c(DV$y[1], TVMTT, ID50P, TVKDEP, TVKDEC, TVINCTRANS, TVINCPROL)))
			fl <- newton.results$p
			MAP.log <- fl[length(fl[,1]),]
			parameters_indiv <- exp(MAP.log) #converged parameter
			parameters_indiv <- c(BASE = parameters_indiv[1],
														MTT = parameters_indiv[2],
														POWER = TVPOWER,
														ID50P = parameters_indiv[3],
														KDEP = parameters_indiv[4],
														KDEC = parameters_indiv[5],
														INCTRANS = parameters_indiv[6],
														INCPROL = parameters_indiv[7])
		}else{
			newton.results <- newton(ll2, x0=log(c(DV$y[1], TVMTT, ID50P, TVTVKDEP, TVKDEC)))
			fl <- newton.results$p
			MAP.log <- fl[length(fl[,1]),]
			parameters_indiv <- exp(MAP.log) #converged parameter
			parameters_indiv <- c(BASE = parameters_indiv[1],
														MTT = parameters_indiv[2],
														POWER = TVPOWER,
														ID50P = parameters_indiv[3],
														KDEP = parameters_indiv[4],
														KDEC = parameters_indiv[5],
														INCTRANS = TVINCTRANS,
														INCPROL = TVINCPROL)
		}
	
		
		
		
		##simulation based on individualized parameters
		
		state <- c(A_1 = 0, 
							 A_2 = 0, 
							 A_3 = parameters_indiv[["BASE"]],
							 A_4 = parameters_indiv[["BASE"]],
							 A_5 = parameters_indiv[["BASE"]],
							 A_6 = parameters_indiv[["BASE"]],
							 A_7 = parameters_indiv[["BASE"]])
		
		out.indiv <- sim(state,totaldose,parameters_indiv,maxtime)
		out <- cbind(out, indiv = out.indiv[,"A_7"])
		
		
		### MCMC
		G = newton.results$G
		H = newton.results$H
		COV.map <- solve(H) 
		
		mcmc.smpl <- MCMC(function(x) -ll2(x), 1000, MAP.log, scale = COV.map/1.5,
				 adapt = F, acc.rate = 0.25)
		mat.resampled <- unique(mcmc.smpl$samples)
		
		### SIR
		if(0){
		# proposal sampling NAP
		
		G = newton.results$G
		H = newton.results$H
		COV.map <- solve(H) 
		
		Nrand <- 100
		set.seed(1111)
		
		mat.resampled <- matrix(NA, nrow = 0, ncol = 7)
		#repeat 
		count <- 0
		while(nrow(mat.resampled) < 100){
			
			count <- count + 1
			print(count)
			proposal.prm.log <- rmvnorm(Nrand,mean =  MAP.log , sigma = COV.map)
			proposal.prob <- dmvnorm(x = proposal.prm.log, mean =  MAP.log , sigma = COV.map )
			
			weight.unnorm <- apply(proposal.prm.log, 1, ll2)
			
			weight.unnorm <- exp(-weight.unnorm)
			weight.unnorm <- weight.unnorm/proposal.prob
			
			weight.norm <- weight.unnorm/sum(weight.unnorm)
			
			idx.resample <- runif(100, 0,1) <= weight.norm
			
			mat.resampled <- rbind(mat.resampled, proposal.prm.log[idx.resample,])
		}
		
		
		}
		
		
		
		for(i in 1:nrow(mat.resampled)){
			parameters_resample <- exp(mat.resampled[i,])
			if(GCSFIIV){
				parameters_resample<- c(BASE = parameters_resample[1],
																MTT = parameters_resample[2],
																POWER = TVPOWER,
																ID50P = parameters_resample[3],
																KDEP = parameters_resample[4],
																KDEC = parameters_resample[5],
																INCTRANS = parameters_resample[6],
																INCPROL = parameters_resample[7])
			}else{
				parameters_resample<- c(BASE = parameters_resample[1],
																MTT = parameters_resample[2],
																POWER = TVPOWER,
																ID50P = parameters_resample[3],
																KDEP = parameters_resample[4],
																KDEC = parameters_resample[5],
																INCTRANS = TVINCTRANS,
																INCPROL = TVINCPROL)
			}
		
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = parameters_resample[["BASE"]],
								 A_4 = parameters_resample[["BASE"]],
								 A_5 = parameters_resample[["BASE"]],
								 A_6 = parameters_resample[["BASE"]],
								 A_7 = parameters_resample[["BASE"]])
			
			out.resample <- sim(state,totaldose, parameters_resample,maxtime)
			out <- cbind(out, indiv.resmpl = out.resample[,"A_7"])
			colnames(out)[ncol(out)] <- paste0("indiv.resmpl",i)
			
		}
		
		
		
		
		
		##simulation based on individualized parameters
		
		
		### MCMC
		
		# pro
		
		
		
		if(0){ 
			Nrand <- 1000
			proposal.eta <- data.frame(ETA_BASE = rep(NA,Nrand),
																 ETA_MTT = rep(NA,Nrand),
																 ETA_ID50 = rep(NA,Nrand),
																 ETA_KDEP = rep(NA,Nrand),
																 ETA_KDEC = rep(NA,Nrand),
																 ETA_TRANS = rep(NA,Nrand),
																 ETA_PROL = rep(NA,Nrand),
																 prob = rep(NA,Nrand))
			
			proposal.prob <- data.frame(ETA_BASE = rep(NA,Nrand),
																	ETA_MTT = rep(NA,Nrand),
																	ETA_ID50 = rep(NA,Nrand),
																	ETA_KDEP = rep(NA,Nrand),
																	ETA_KDEC = rep(NA,Nrand),
																	ETA_TRANS = rep(NA,Nrand),
																	ETA_PROL = rep(NA,Nrand))
			
			set.seed(1111)
			for(i in 1:7){  
				proposal.eta[,i] <- rnorm(Nrand,0, omega[i])
			}
			
			
			for(i in 1:7){
				proposal.prob[,i] <- dnorm(x = proposal.eta[,i],sd = omega[i],mean = 0 )
			}
			
			proposal.eta$prob <- exp(rowSums(log(proposal.prob)))
			proposal.prm <- proposal.eta[,-8]
			
			
			
			for(i in 1:7){  
				proposal.prm[,i ] <- parameters[-3][i]*exp(proposal.prm[,i])
			}
			
			proposal.prm <- as.matrix(proposal.prm)
			colnames(proposal.prm) <- NULL
			weight.unnorm <- apply(log(proposal.prm), 1, ll2)
			
			weight.unnorm <- exp(-weight.unnorm)
			
			
			
			ll2(log(proposal.prm)[1,])
			weight.unnorm <- weight.unnorm
			weight.unnorm <- weight.unnorm/proposal.eta$prob
			weight.norm <- weight.unnorm/sum(weight.unnorm)
			
			
			
			randunif <- runif(100000, 0,1)
			
			idx.resample <- randunif <= weight.norm
		}
		
		
		
		
		
		
	}
	
	
	
	return (out)
	#write(parameters, 'finalpara')
	write.csv(out, 'simulation.csv')
}

runmain_typical <- function(totaldose,  SEX, DM,  maxtime, BASE, DV = NULL, plot.variability ) {
	
	
	GCSFIIV = 0
	
	if(GCSFIIV){
		TVBASE <- 5.324
		BASE <- BASE
		TVMTT <- 4.668
		TVPOWER  <-  0.1953
		TVID50P <-  88.48
		TVKDEP <- 0.03704
		TVKDEC <- 0.1718
		TVINCTRANS <- 3.073
		TVINCPROL <- 0.208
		ID50P <- TVID50P*(1+0.4596*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2374,  MTT=0.1608, ID50P =0.2722, KDEP=0.8246,KDEC= 0.9515,INCTRANS= 0.3388,INCPROL= 0.3388)
		sigadd <- 0.529
		sigprop <- 0.2941
		
	}else{
		TVBASE <- 5.341
		BASE <- BASE
		TVMTT <- 4.639
		TVPOWER  <-  0.1882
		TVID50P <-  88.88
		TVKDEP <- 0.03261
		TVKDEC <- 0.1877
		TVINCTRANS <- 3.502
		TVINCPROL <- 0.2172
		ID50P <- TVID50P*(1+0.4849*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2301, MTT= 0.1616,ID50P = 0.2332,KDEP= 0.9557,  KDEC= 0.92)
		sigadd <- 0.5271
		sigprop <- 0.2998
		
	}
	
	
	######################################
	#####population mean simulation ######
	######################################
	
	parameters <- c(BASE = BASE,
									MTT = TVMTT,
									POWER = TVPOWER,
									ID50P = ID50P,
									KDEP = TVKDEP,
									KDEC = TVKDEC,
									INCTRANS = TVINCTRANS,
									INCPROL =TVINCPROL)
	
	state <- c(A_1 = 0, 
						 A_2 = 0, 
						 A_3 = BASE,
						 A_4 = BASE,
						 A_5 = BASE,
						 A_6 = BASE,
						 A_7 = BASE)
	
	
	out <- sim(state,totaldose, parameters, maxtime)
	colnames(out)[8] <- "typical"
	
	
	#############################
	###### variabilbity #########
	#############################
	
	if(as.numeric(plot.variability) == 1){
		
		Nrand <- 500
		set.seed(1111)
		
		if(GCSFIIV){
			prm.eta <- data.frame(ETA_BASE = rep(NA,Nrand),
														ETA_MTT = rep(NA,Nrand),
														ETA_ID50 = rep(NA,Nrand),
														ETA_KDEP = rep(NA,Nrand),
														ETA_KDEC = rep(NA,Nrand),
														ETA_TRANS = rep(NA,Nrand),
														ETA_PROL = rep(NA,Nrand))
			for(i in 1:7){  
				
				prm.eta[,i ] <- parameters[-3][i]*exp(rnorm(Nrand,0, omega[[i]]))
			}
			
		}else{
			prm.eta <- data.frame(ETA_BASE = rep(NA,Nrand),
														ETA_MTT = rep(NA,Nrand),
														ETA_ID50 = rep(NA,Nrand),
														ETA_KDEP = rep(NA,Nrand),
														ETA_KDEC = rep(NA,Nrand))
			for(i in 1:5){  
				prm.eta[,i ] <- parameters[-3][i]*exp(rnorm(Nrand,0, omega[[i]]))
			}
		}
		
		
		###simulation for interindividual variability
		for(idx.var in 1:Nrand){
			if(GCSFIIV){
				parameters.eta <- c(BASE = BASE,
														MTT = prm.eta[idx.var,2],
														POWER = TVPOWER,
														ID50P = prm.eta[idx.var,3],
														KDEP = prm.eta[idx.var,4],
														KDEC = prm.eta[idx.var,5],
														INCTRANS = prm.eta[idx.var,6],
														INCPROL =prm.eta[idx.var,7])
			}else{
				parameters.eta <- c(BASE = BASE,
														MTT = prm.eta[idx.var,2],
														POWER = TVPOWER,
														ID50P = prm.eta[idx.var,3],
														KDEP = prm.eta[idx.var,4],
														KDEC = prm.eta[idx.var,5],
														INCTRANS = TVINCTRANS,
														INCPROL =TVINCPROL)
			}
			
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = BASE,
								 A_4 = BASE,
								 A_5 = BASE,
								 A_6 = BASE,
								 A_7 = BASE)
			
			
			
			temp.out <- sim(state,totaldose, parameters.eta, maxtime)
			out <- cbind(out, var = temp.out[,"A_7"])
		}
		colnames(out)[9:508] <- paste0("var", 1:500)
	}
	
	return(out)	
}


##individual
runmain_indiv_estim <- function(totaldose,  SEX, DM,  maxtime, BASE, DV = NULL, estim.variability ) {
	
	
	GCSFIIV = 0
	BASE.ESTIM = 0
	
	if(GCSFIIV){
		TVBASE <- 5.324
		BASE <- BASE
		TVMTT <- 4.668
		TVPOWER  <-  0.1953
		TVID50P <-  88.48
		TVKDEP <- 0.03704
		TVKDEC <- 0.1718
		TVINCTRANS <- 3.073
		TVINCPROL <- 0.208
		ID50P <- TVID50P*(1+0.4596*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2374,  MTT=0.1608, ID50P =0.2722, KDEP=0.8246,KDEC= 0.9515,INCTRANS= 0.3388,INCPROL= 0.3388)
		sigadd <- 0.529
		sigprop <- 0.2941
		
	}else{
		TVBASE <- 5.341
		BASE <- BASE
		TVMTT <- 4.639
		TVPOWER  <-  0.1882
		TVID50P <-  88.88
		TVKDEP <- 0.03261
		TVKDEC <- 0.1877
		TVINCTRANS <- 3.502
		TVINCPROL <- 0.2172
		ID50P <- TVID50P*(1+0.4849*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2301, MTT= 0.1616,ID50P = 0.2332,KDEP= 0.9557,  KDEC= 0.92)
		sigadd <- 0.5271
		sigprop <- 0.2998
		
	}
	
	
	
	
	###############################################
	#######individualized parameter estimation#####
	###############################################
	
	if(!is.null(DV)){
		
		parameters <- c(BASE = TVBASE,
										MTT = TVMTT,
										POWER = TVPOWER,
										ID50P = ID50P,
										KDEP = TVKDEP,
										KDEC = TVKDEC,
										INCTRANS = TVINCTRANS,
										INCPROL =TVINCPROL)
		
		# Negative Log-likelihood
		llmod <- function(state,totaldose,parameters, maxtime){
			ind_pred <- sim(state,totaldose,parameters, maxtime)
			ind_pred <- ind_pred %>% filter(time %in% DV$time) %>% select(time, A_7) %>% mutate(y = A_7) %>% select(-A_7) %>% unique()
			like <- -dnorm((DV$y-ind_pred$y)/sqrt(sigprop^2*(ind_pred$y)^2+sigadd^2), 0, 1, log=T)
			#SQRT(THETA(9)**2*IPRED**2+THETA(10)**2)
			
			return(like)
		}
		
		ll <- function(parameters) {
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = parameters[["BASE"]],
								 A_4 = parameters[["BASE"]],
								 A_5 = parameters[["BASE"]],
								 A_6 = parameters[["BASE"]],
								 A_7 = parameters[["BASE"]])
			if(GCSFIIV){
				temp.sum <- sum(llmod(state,totaldose,parameters, max(DV$time))) - 
					dnorm(log(parameters["BASE"]), log(TVBASE), omega["BASE"], log=T) -
					dnorm(log(parameters["MTT"]), log(TVMTT), omega["MTT"], log=T) -
					dnorm(log(parameters["ID50P"]), log(ID50P), omega["ID50P"], log=T) -
					dnorm(log(parameters["KDEP"]), log(TVKDEP), omega["KDEP"], log=T) -
					dnorm(log(parameters["KDEC"]), log(TVKDEC), omega["KDEC"], log=T) -
					dnorm(log(parameters["INCTRANS"]), log(TVINCTRANS), omega["INCTRANS"], log=T) -
					dnorm(log(parameters["INCPROL"]), log(TVINCPROL), omega["INCPROL"], log=T)
			}else{
				temp.sum <- sum(llmod(state,totaldose,parameters, max(DV$time))) - 
					dnorm(log(parameters["BASE"]), log(TVBASE), omega["BASE"], log=T) -
					dnorm(log(parameters["MTT"]), log(TVMTT), omega["MTT"], log=T) -
					dnorm(log(parameters["ID50P"]), log(ID50P), omega["ID50P"], log=T) -
					dnorm(log(parameters["KDEP"]), log(TVKDEP), omega["KDEP"], log=T) -
					dnorm(log(parameters["KDEC"]), log(TVKDEC), omega["KDEC"], log=T) 
			}
			
		
		
			return(temp.sum)
		}
		
		ll2 <- function(x){
			if(GCSFIIV){
				if(BASE.ESTIM){
					parameters <- c(BASE = exp(x[1]),
													MTT = exp(x[2]),
													POWER = TVPOWER,
													ID50P = exp(x[3]),
													KDEP = exp(x[4]),
													KDEC = exp(x[5]),
													INCTRANS = exp(x[6]),
													INCPROL = exp(x[7]))
				}else{
					parameters <- c(BASE = BASE,
													MTT = exp(x[1]),
													POWER = TVPOWER,
													ID50P = exp(x[2]),
													KDEP = exp(x[3]),
													KDEC = exp(x[4]),
													INCTRANS = exp(x[5]),
													INCPROL = exp(x[6]))
				}
			}else{
				if(BASE.ESTIM){
					parameters <- c(BASE = exp(x[1]),
													MTT = exp(x[2]),
													POWER = TVPOWER,
													ID50P = exp(x[3]),
													KDEP = exp(x[4]),
													KDEC = exp(x[5]),
													INCTRANS = TVINCTRANS,
													INCPROL = TVINCPROL)
				}else{
					parameters <- c(BASE = BASE,
													MTT = exp(x[1]),
													POWER = TVPOWER,
													ID50P = exp(x[2]),
													KDEP = exp(x[3]),
													KDEC = exp(x[4]),
													INCTRANS = TVINCTRANS,
													INCPROL = TVINCPROL)
				}
			
			}
			
			return(ll(parameters))
		}
		#ll2(log(c(5.324, 4.668, 88.48, 0.03704, 0.1718, 3.073, 0.208)))
		
		
		#newton.results <- newton(ll2, x0=log(c(5.324, 4.668, ID50P, 0.03704, 0.1718, 3.073, 0.208)))
		if(GCSFIIV){
			if(BASE.ESTIM){
				newton.results <- newton(ll2, x0=log(c(DV$y[1], TVMTT, ID50P, TVKDEP, TVKDEC, TVINCTRANS, TVINCPROL)))
				fl <- newton.results$p
				MAP.log <- fl[length(fl[,1]),]
				parameters_indiv <- exp(MAP.log) #converged parameter
				parameters_indiv <- c(BASE = parameters_indiv[1],
															MTT = parameters_indiv[2],
															POWER = TVPOWER,
															ID50P = parameters_indiv[3],
															KDEP = parameters_indiv[4],
															KDEC = parameters_indiv[5],
															INCTRANS = parameters_indiv[6],
															INCPROL = parameters_indiv[7])
			}else{
				newton.results <- newton(ll2, x0=log(c(TVMTT, ID50P, TVKDEP, TVKDEC, TVINCTRANS, TVINCPROL)))
				fl <- newton.results$p
				MAP.log <- fl[length(fl[,1]),]
				parameters_indiv <- exp(MAP.log) #converged parameter
				parameters_indiv <- c(BASE = BASE,
															MTT = parameters_indiv[1],
															POWER = TVPOWER,
															ID50P = parameters_indiv[2],
															KDEP = parameters_indiv[3],
															KDEC = parameters_indiv[4],
															INCTRANS = parameters_indiv[5],
															INCPROL = parameters_indiv[6])
			}
		}else{
			if(BASE.ESTIM){
				newton.results <- newton(ll2, x0=log(c(DV$y[1], TVMTT, ID50P, TVKDEP, TVKDEC)))
				fl <- newton.results$p
				MAP.log <- fl[length(fl[,1]),]
				parameters_indiv <- exp(MAP.log) #converged parameter
				parameters_indiv <- c(BASE = parameters_indiv[1],
															MTT = parameters_indiv[2],
															POWER = TVPOWER,
															ID50P = parameters_indiv[3],
															KDEP = parameters_indiv[4],
															KDEC = parameters_indiv[5],
															INCTRANS = TVINCTRANS,
															INCPROL = TVINCPROL)
			}else{
				newton.results <- newton(ll2, x0=log(c( TVMTT, ID50P, TVKDEP, TVKDEC)))
				fl <- newton.results$p
				MAP.log <- fl[length(fl[,1]),]
				parameters_indiv <- exp(MAP.log) #converged parameter
				parameters_indiv <- c(BASE = BASE,
															MTT = parameters_indiv[1],
															POWER = TVPOWER,
															ID50P = parameters_indiv[2],
															KDEP = parameters_indiv[3],
															KDEC = parameters_indiv[4],
															INCTRANS = TVINCTRANS,
															INCPROL = TVINCPROL)
			}
			
		}
		
		
		
		##simulation based on individualized parameters
		if(as.numeric(estim.variability) == 1){
	
			### MCMC
			G = newton.results$G
			H = newton.results$H
			COV.map <- solve(H) 
			COV.map <- (COV.map + t(COV.map))/2
			
			COV.map[diag(COV.map)<0,diag(COV.map)<0] <- abs(COV.map[diag(COV.map)<0,diag(COV.map)<0])
			
			if(0){
			mcmc.smpl <- MCMC(function(x) -ll2(x), 1000, MAP.log, scale = COV.map/1.5,
												adapt = F, acc.rate = 0.25)
			mcmc.smpl <-MCMC(function(x) -ll2(x), 1000, MAP.log, scale = abs(diag(COV.map)),
											 adapt = F, acc.rate = 0.25)
			

			}
			
			mcmc.smpl <- MCMC(function(x) -ll2(x), 1000, MAP.log, scale = COV.map,
												adapt = F, acc.rate = 0.25)
			mat.resampled <- unique(mcmc.smpl$samples)
		}
		
	
	}
		
	
		
	
	if(as.numeric(estim.variability) == 1){
		return(list(parameters_indiv = parameters_indiv,
								mat.resampled = mat.resampled))
	}else{
		return(list(parameters_indiv = parameters_indiv))
	}
	
}


runmain_indiv_pred <- function(totaldose,  SEX, DM,  maxtime, BASE, DV = NULL, list.estimated.prms, plot.variability ) {
	
	
	GCSFIIV = 0
	BASE.ESTIM = 0
	
	if(GCSFIIV){
		TVBASE <- 5.324
		BASE <- BASE
		TVMTT <- 4.668
		TVPOWER  <-  0.1953
		TVID50P <-  88.48
		TVKDEP <- 0.03704
		TVKDEC <- 0.1718
		TVINCTRANS <- 3.073
		TVINCPROL <- 0.208
		ID50P <- TVID50P*(1+0.4596*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2374,  MTT=0.1608, ID50P =0.2722, KDEP=0.8246,KDEC= 0.9515,INCTRANS= 0.3388,INCPROL= 0.3388)
		sigadd <- 0.529
		sigprop <- 0.2941
		
	}else{
		TVBASE <- 5.341
		BASE <- BASE
		TVMTT <- 4.639
		TVPOWER  <-  0.1882
		TVID50P <-  88.88
		TVKDEP <- 0.03261
		TVKDEC <- 0.1877
		TVINCTRANS <- 3.502
		TVINCPROL <- 0.2172
		ID50P <- TVID50P*(1+0.4849*as.numeric(DM))*(1-0.3338*as.numeric(SEX))
		omega <- c(BASE=0.2301, MTT= 0.1616,ID50P = 0.2332,KDEP= 0.9557,  KDEC= 0.92)
		sigadd <- 0.5271
		sigprop <- 0.2998
		
	}
	
	

	######################################
	#####MAP simulation ######
	######################################
	
	parameters <- c(BASE = BASE,
									MTT = TVMTT,
									POWER = TVPOWER,
									ID50P = ID50P,
									KDEP = TVKDEP,
									KDEC = TVKDEC,
									INCTRANS = TVINCTRANS,
									INCPROL =TVINCPROL)
	
	
	###############################################
	#######individualized parameter estimation#####
	###############################################
	
		
	parameters_indiv <- list.estimated.prms$parameters_indiv

	
		parameters_indiv <- c(BASE = parameters_indiv[[1]],
													MTT = parameters_indiv[[2]],
													POWER = parameters_indiv[[3]],
													ID50P = parameters_indiv[[4]],
													KDEP = parameters_indiv[[5]],
													KDEC = parameters_indiv[[6]],
													INCTRANS = parameters_indiv[[7]],
													INCPROL = parameters_indiv[[8]])
		
	state <- c(A_1 = 0, 
						 A_2 = 0, 
						 A_3 = parameters_indiv[["BASE"]],
						 A_4 = parameters_indiv[["BASE"]],
						 A_5 = parameters_indiv[["BASE"]],
						 A_6 = parameters_indiv[["BASE"]],
						 A_7 = parameters_indiv[["BASE"]])
	
	out.indiv <- sim(state,totaldose,parameters_indiv,maxtime)
	colnames(out.indiv)[8] <- "indiv"
	#out <- cbind(out, indiv = out.indiv[,"A_7"])
	
	
		##simulation based on individualized parameters
	if( as.numeric(plot.variability) == 1){
		mat.resampled <- list.estimated.prms$mat.resampled
		for(i in 1:nrow(mat.resampled)){
			parameters_resample <- exp(mat.resampled[i,])
			if(GCSFIIV){
				if(BASE.ESTIM){
					parameters_resample<- c(BASE = parameters_resample[1],
																	MTT = parameters_resample[2],
																	POWER = TVPOWER,
																	ID50P = parameters_resample[3],
																	KDEP = parameters_resample[4],
																	KDEC = parameters_resample[5],
																	INCTRANS = parameters_resample[6],
																	INCPROL = parameters_resample[7])
				}else{
					parameters_resample<- c(BASE = BASE,
																	MTT = parameters_resample[1],
																	POWER = TVPOWER,
																	ID50P = parameters_resample[2],
																	KDEP = parameters_resample[3],
																	KDEC = parameters_resample[4],
																	INCTRANS = parameters_resample[5],
																	INCPROL = parameters_resample[6])
				}
			}else{
				if(BASE.ESTIM){
					parameters_resample<- c(BASE = parameters_resample[1],
																	MTT = parameters_resample[2],
																	POWER = TVPOWER,
																	ID50P = parameters_resample[3],
																	KDEP = parameters_resample[4],
																	KDEC = parameters_resample[5],
																	INCTRANS = TVINCTRANS,
																	INCPROL = TVINCPROL)	
				}else{
					parameters_resample<- c(BASE = BASE,
																	MTT = parameters_resample[1],
																	POWER = TVPOWER,
																	ID50P = parameters_resample[2],
																	KDEP = parameters_resample[3],
																	KDEC = parameters_resample[4],
																	INCTRANS = TVINCTRANS,
																	INCPROL = TVINCPROL)
				}
				
			}
			
			state <- c(A_1 = 0, 
								 A_2 = 0, 
								 A_3 = parameters_resample[["BASE"]],
								 A_4 = parameters_resample[["BASE"]],
								 A_5 = parameters_resample[["BASE"]],
								 A_6 = parameters_resample[["BASE"]],
								 A_7 = parameters_resample[["BASE"]])
			
			out.resample <- sim(state,totaldose, parameters_resample,maxtime)
			out.indiv <- cbind(out.indiv, indiv.resmpl = out.resample[,"A_7"])
			colnames(out.indiv)[ncol(out.indiv)] <- paste0("indiv.resmpl",i)
		}
	}
		
	
	
	
	return (out.indiv)
}


