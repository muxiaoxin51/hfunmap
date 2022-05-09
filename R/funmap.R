#' hfun_load_data    Initial database load
#'
#' @param address.p charactor like "/home/C/data/pheno.csv"
#' @param address.g charactor like "/home/C/data/SNP.csv",missing value is numeric -9
#' @param Times numeric of times like c(1,2,3,6,14,19,25,48)
#' @param ranks numerc the number of traits
#' @param missing logical whether geno database have missing value, default is False
#' @param log.p logical whether transform pheno to log10.pheno, default is True
#' @param fit.check logical whether do individual curve fit; return mean R^2
#'
#' @return list dat.p,dat.g,dat.t,dat.r,dat.m
#' @export
hfun_load_data <- function(address.p,address.g,Times,ranks=NULL,missing=FALSE,
                      log.p=TRUE,fit.check=FALSE){
  n.times <- length(Times)
  input.p <- utils::read.csv(address.p,row.names = 1)
  input.g <- utils::read.csv(address.g,row.names = 1)
  if(is.null(ranks)){
    ranks <- dim(input.p)[2]/n.times
  }
  if(log.p==TRUE){
    if(min(log.p)==0){input.p <- input.p + 1}
    input.p <- log10(input.p)
  }
  n.iden <- dim(input.p)[1]
  ng <- length(table(as.matrix(input.g)))
  if(missing==TRUE){
    ng <- ng - 1
  }
  if(fit.check==FALSE){
    return(list(dat.p=input.p,dat.g=input.g,dat.t=Times,dat.r=ranks,
                dat.m=missing,dat.ng=ng))
  }else{
    test <- t(apply(as.matrix(input.p),1,function(c){hfun_iden_fit(x=c,Times=Times,ranks=ranks)}))
    mean.fit <- colMeans(test)[seq(4,dim(test)[2],4)]
    cat(mean.fit,"\n")
    return(list(dat.p=input.p,dat.g=input.g,dat.t=Times,dat.r=ranks,
                dat.m=missing,dat.ng=ng))
  }
}



#' hfun_iden_fit individual curve fit
#'
#' @param x vector of individual pheno
#' @param Times numeric of times like c(1,2,3,6,14,19,25,48)
#' @param ranks numerc the number of traits
#'
#' @return vector, which is composed of ranks part of curve_par a,b,r and R^2
#' @export

hfun_iden_fit <- function(x,Times,ranks){
  res <- c()
  ols_fit <- function(par,real,Times){
    exp <- hfun_g_t(par,Times)
    ols <- sum((exp - real)^2)
    return(ols)
  }
  for(i in 1:ranks){
    input <- x[seq(i,length(x),ranks)]
    tpp <- DEoptim::DEoptim(fn=ols_fit,lower =rep(1e-3,3) ,upper =rep(10*max(input),3) ,
                            control = list(trace=F),real=input,Times=Times)
    tpp <- stats::optim(tpp$optim$bestmem,ols_fit,real=input,Times=Times,method = "BFGS",control = list(maxit=5000))
    R2 <- 1 - (sum((hfun_g_t(tpp$par,Times)-input)^2)/sum((x - mean(x))^2))
    res <- c(res,tpp$par,R2)
  }
  return(res)
}

#' hfun_get_ini Initialization of fitting mean and covariance matrix parameters
#'
#' @param phe matrix:input phenotype data
#' @param Times numeric of times like c(1,2,3,6,14,19,25,48)
#' @param ranks numerc the number of traits
#' @param control list include
#' cur_lower        numeric of curve parameters lower set, default is rep(1e-3,3)
#' cur_upper        numeric of curve parameters upper set, default is rep(1000,3)
#' cur_control      list of curve DEoptim control, default is list(itermax=500,trace=F), more details see ?DEoptim
#' cov_lower        numeric of cov-matrix parameters lower set, default is 1e-3
#' cov_upper        numeric of cov-matrix parameters upper set, default is 100
#' cov_control      list of cov-matrix DEoptim control, default is list(itermax=500,trace=F), more details see ?DEoptim
#'
#' @return numeric parameters of curve and cov-matrix
#' @export
hfun_get_ini <- function(phe,Times,ranks,control=list(cur_lower=rep(1e-3,3),cur_upper=rep(1000,3),cur_control=list(itermax=500,trace=F),
                         cov_lower=1e-3,cov_upper=100,cov_control=list(itermax=500,trace=F))){
  curve_fit <- function(par,y,Times){
    n <- dim(y)[1]
    miu <- hfun_g_t(par,Times)
    y.delt <- y - matrix(rep(miu,n),nrow = n,byrow = TRUE)
    res <- sum(y.delt^2)/n
    return(res)
  }
  tpp <- list()
  tppp <- list()
  for(i in 1:ranks){
    tpp[[i]] <- DEoptim::DEoptim(fn=curve_fit,lower = control$cur_lower,upper = control$cur_upper,control = control$cur_control,y = phe[,seq(i,dim(phe)[2],ranks)],Times=Times)
    tppp[[i]] <- try(stats::optim(tpp[[i]]$optim$bestmem,curve_fit,y=phe[,seq(i,dim(phe)[2],ranks)],Times=Times,method = "BFGS",
                           control = list(maxit=5000,trace=F)),silent = T)
  }
  #tpp <- DEoptim(fn=curve_fit,lower = rep(control$cur_lower,ranks),upper = rep(control$cur_upper,ranks),control = control$cur_control,y = phe,Times=Times)
  #tppp <- try(optim(tpp$optim$bestmem,curve_fit,y=phe,Times=Times,method = "BFGS",
  #                  control = list(maxit=5000,trace=F)),silent = T)
  if(length(which(c(unlist(lapply(tppp,function(c){class(c)})))=="try-error"))>0){
    cat("There is something wrong in curve fit parameter initial!\n")
    return("error")
  }else{
    curve_par <- as.numeric(unlist(lapply(tppp,function(x){x$par})))
    y.resd <- phe - matrix(rep(hfun_get_mu(curve_par,Times),dim(phe)[1]),nrow = dim(phe)[1],byrow = TRUE)
    cov_fit <- function(par,resd,Times){
      m <- length(Times)
      cov_mat <- hfun_get_invsigma(par,Times,dim(resd)[2]/m)
      pv <- log(1/(2*pi)^(m/2)*sqrt(abs(cov_mat$det))) + ((-1/2)*rowSums(resd %*% cov_mat$sigma_1 * resd))
      Lh <- sum(pv)/dim(resd)[1]
      return(-Lh)
    }
    n_par_cov <- (3 + ranks)*ranks/2
    cpp <- list()
    cppp <- list()
    for(i in 1:ranks){
      cpp[[i]] <- DEoptim::DEoptim(fn=cov_fit,lower = rep(1e-8,2),upper = c(1,1000),control = control$cov_control,
                     resd=y.resd[,seq(i,dim(y.resd)[2],ranks)],Times=Times)
      #cppp_lower <- cpp[[i]]$optim$bestmem * 0.95
      #cppp_upper <- cpp[[i]]$optim$bestmem * 1.05
      #cppp[[i]] <- try(optim(cpp[[i]]$optim$bestmem,cov_fit,resd=y.resd[,seq(i,dim(y.resd)[2],ranks)],Times=Times,method = "L-BFGS-B",
      #                       lower = cppp_lower,upper = cppp_upper,control = list(maxit=5000,trace=F)))
      cppp[[i]] <- try(stats::optim(cpp[[i]]$optim$bestmem,cov_fit,resd=y.resd[,seq(i,dim(y.resd)[2],ranks)],Times=Times,method = "BFGS",
                             control = list(maxit=5000,trace=F)))
    }
    cpp_par <- as.numeric(matrix(unlist(lapply(cppp,function(x){x$par})),ncol =  2,byrow = T))
    cpp_par_lower <- c(cpp_par*0.95,rep(1e-8,ranks*(ranks-1)/2))
    cpp_par_upper <- c(cpp_par*1.05,rep(1,ranks*(ranks-1)/2))

    #cppn <- DEoptim(fn=cov_fit,lower = cpp_par_lower,upper = cpp_par_upper,control = control$cov_control,
    #               resd=y.resd,Times=Times)
    iter <- 0
    cpppn <- list(par=cpp_par,value=-1)
    while(cpppn$value<0 && iter <6){
      cpppn <- try(stats::optim(c(cpp_par,rep(0.5,ranks*(ranks-1)/2)),cov_fit,resd=y.resd,Times=Times,method = "BFGS",
                         control = list(maxit=5000,trace=F,parscale=rep(.1^iter,n_par_cov),abstol=TRUE)))
      iter <- iter + 1
      if(class(cpppn)=="try-error"){cpppn <- list(par=cpp_par,value=-1)}
    }
    iter <- iter - 1
    if(class(cpppn)=="try-error"){
      cat("There is something wrong in cov fit parameter initial!\n")
      return("error")
    }else{
      all_par <- c(curve_par,cpppn$par)
      return(all_par)
    }
  }
}

#' hfun_lr_cal LR calculate function
#'
#' @param format_data data from hfun_load_data
#' @param interval numeric the serial number of genotype you want to calculate eg.c(1:10) or c(1,5,9,12,56)
#' @param useinv  logical whether use inverse matrix, default is TRUE
#' @param usecpp  logical whether use Rcpp code, default is TRUE
#' @param ini_control list,default is NULL, include
#' cur_lower        numeric of curve parameters lower set, default is rep(1e-3,3)
#' cur_upper        numeric of curve parameters upper set, default is rep(1000,3)
#' cur_control      list of curve DEoptim control, default is list(itermax=500,trace=F), more details see ?DEoptim
#' cov_lower        numeric of cov-matrix parameters lower set, default is 1e-3
#' cov_upper        numeric of cov-matrix parameters upper set, default is 100
#' cov_control      list of cov-matrix DEoptim control, default is list(itermax=500,trace=F), more details see ?DEoptim
#' @param LR_only logical whether output only include LR value, default is FALSE
#'
#' @return matrix of lr_calculate form example c(LR,H1_mle$value,H1_mle$par,H0_mle$value,H0_mle$par)
#' @export
hfun_lr_cal <- function(format_data,interval,useinv=TRUE,usecpp=TRUE,ini_control=NULL,LR_only=FALSE){
  pheno <- as.matrix(format_data$dat.p)
  geno <- as.matrix(format_data$dat.g)
  times <- format_data$dat.t
  ranks <- format_data$dat.r
  dat_missing <- format_data$dat.m
  if(is.null(ini_control)){ini_control <- list(cur_lower=rep(1e-3,3),cur_upper=rep(1000,3),cur_control=list(itermax=500,trace=F),
                                            cov_lower=1e-3,cov_upper=100,cov_control=list(itermax=500,trace=F))}
  interval <- sort(interval)
  ng <- format_data$dat.ng
  res_col <- (3*ranks+(3+ranks)*ranks/2) + (3*ranks*ng+(3+ranks)*ranks/2) + 3
  res <- matrix(data = 0,nrow = max(interval),ncol = res_col)

  if(dat_missing==TRUE){
    for(i in interval){
      nSNP <- geno[i,]
      missings <- which(nSNP==-9)
      if(length(missings)>0){
        SNP1 <- nSNP[-(missings)]
        pheno1 <- pheno[-(missings),]
      }else{
        SNP1 <- nSNP
        pheno1 <- pheno
      }
      nphe <- pheno1
      NSNP <- as.numeric(SNP1)
      iterm <- 0
      while(iterm<10){
        ini_par <- try(hfun_get_ini(phe = nphe,Times = times,ranks = ranks,control = ini_control),silent = T)
        if(class(ini_par)=="try-error"){
          iterm <- iterm + 1
          next
        }else if(min(ini_par)<0){
          iterm <- iterm + 1
          next
        }else{
          break
        }
      }
      if(class(ini_par)=="try-error"){
        res[i,1] <- NA
        cat("SNP",i,": LR=",NA,"para =",NA,"\n")
        next
      }else if(min(ini_par)<0){
        res[i,1] <- NA
        cat("SNP",i,": LR=",NA,"para =",NA,"\n")
        next
      }
      H0_lower <- ini_par * 0.95
      H0_upper <- ini_par * 1.05
      H0_res <- try(stats::optim(par = ini_par,fn=hfun_H0_mle,N_phe=nphe,pheT=times,ranks=ranks,
                      method = "L-BFGS-B",lower = H0_lower,upper = H0_upper,
                      control = list(maxit=5000,trace=F,REPORT=1,parscale=rep(.1^0,length(ini_par)))),silent = TRUE)
      if(class(H0_res)=="try-error"){
        res[i,1] <- NA
        cat("SNP",i,": LR=",NA,"para =",NA,"\n")
        next
      }
      H0_value <- H0_res$value
      nnd <- length(table(NSNP))
      differ.gene <- names(table(NSNP))
      pheno_list <- list()
      H1_ini_par <- c(rep(H0_res$par[1:(3*ranks)],nnd),H0_res$par[-(1:(3*ranks))])
      H1_lower <- H1_ini_par * 0.99
      H1_upper <- H1_ini_par * 1.01
      for(s in 1:nnd){
        class.n <- which(NSNP==differ.gene[s])
        pheno_list[[s]] <- as.matrix(nphe[class.n,])
      }
      H1_res <- try(stats::optim(par = H1_ini_par,fn=hfun_H1_mle,pheno=pheno_list,Times=times,ranks=ranks,
                          ng=nnd,method = "L-BFGS-B",lower = H1_lower,upper = H1_upper,control = list(maxit=5000,
                          trace=F,REPORT=1)),silent = TRUE)
      if(class(H1_res)=="try-error"){
        res[i,1] <- NA
        cat("SNP",i,": LR=",NA,"para =",NA,"\n")
        next
      }
      LR <- 2 * (H0_value - H1_res$value)
      all_para <- as.numeric(c(LR,H1_res$value,H1_res$par,H0_res$value,H0_res$par))
      res[i,1:length(all_para)] <- all_para
      cat("SNP",i,": LR =",LR,"\n")
    }
  }else{
    ini_par <- try(hfun_get_ini(phe = pheno,Times = times,ranks = ranks,control = ini_control),silent = T)
    if(class(ini_par)=="try-error"){
        cat("something is wrong in HO parameters initialize!\n")
        return("error!\n")
    }
    H0_lower <- ini_par * 0.95
    H0_upper <- ini_par * 1.05
    H0_res <- try(stats::optim(par = ini_par,fn=hfun_H0_mle,N_phe=pheno,pheT=times,ranks=ranks,
                        method = "L-BFGS-B",lower = H0_lower,upper = H0_upper,
                        control = list(maxit=5000,trace=F,REPORT=1,parscale=rep(.1^0,length(ini_par)))),silent = TRUE)
    if(class(H0_res)=="try-error"){
      cat("something is wrong in H0 optim!\n")
      return("error!\n")
    }
    H0_value <- H0_res$value
    for(i in interval){
      nSNP <- geno[i,]
      nphe <- pheno
      NSNP <- as.numeric(nSNP)
      nnd <- length(table(NSNP))
      differ.gene <- names(table(NSNP))
      pheno_list <- list()
      H1_ini_par <- c(rep(H0_res$par[1:(3*ranks)],nnd),H0_res$par[-(1:(3*ranks))])
      H1_lower <- H1_ini_par * 0.99
      H1_upper <- H1_ini_par * 1.05
      for(s in 1:nnd){
        class.n <- which(NSNP==differ.gene[s])
        pheno_list[[s]] <- as.matrix(nphe[class.n,])
      }

      H1_res <- try(stats::optim(par = H1_ini_par,fn=hfun_H1_mle,pheno=pheno_list,Times=times,ranks=ranks,
                          ng=nnd,method = "L-BFGS-B",lower = H1_lower,upper = H1_upper,control = list(maxit=5000,
                                                                                                      trace=F,REPORT=1)),silent = TRUE)
      if(class(H1_res)=="try-error"){
        res[i,1] <- NA
        cat("SNP",i,": LR=",NA,"para =",NA,"\n")
        next
      }
      LR <- 2 * (H0_value - H1_res$value)
      all_para <- as.numeric(c(LR,H1_res$value,H1_res$par,H0_res$value,H0_res$par))
      res[i,1:length(all_para)] <- all_para
      cat("SNP",i,": LR =",LR,"\n")
    }
  }
  output <- res[interval,]
  if(LR_only==TRUE){
    output <- output[,1]
  }
  return(output)
}



#' hfun_permutation do permutation test
#'
#' @param format_data data from hfun_load_data
#' @param permu_times int, times you want to do permutation test,default is 5
#' @param interval numeric the serial number of genotype you want to calculate eg.c(1:10) or c(1,5,9,12,56)
#' @param cutp double, cutlevel you want, default is 0.05(P<=0.05)
#'
#' @return list(cut1,cut2), cut1 is calculate by sort(max(LR for single permutation))cutp;cut2 is calculate by max(sort(LR for single permutation)cutp)
#' @export
hfun_permutation <- function(format_data,permu_times=5,interval,cutp=0.05){
  pheno <- as.matrix(format_data$dat.p)
  geno <- as.matrix(format_data$dat.g)
  times <- format_data$dat.t
  ranks <- format_data$dat.r
  dat_missing <- format_data$dat.m
  sample_num <- dim(pheno)[1]
  c <- matrix(NA,nrow = permu_times,ncol = sample_num)
  for(i in 1:permu_times){
    c[i,]<- sample(c(1:sample_num),sample_num)
  }
  lr_all <- sapply(c(1:permu_times),function(x){
    format_data$dat.p <- pheno[c[i,],];
    y <- hfun_lr_cal(format_data = format_data,interval = interval,LR_only = TRUE)
    return(y)
    })
  nn <- round(permu_times * cutp)
  if(nn<1){nn <- 1}
  cutlevel1 <- sort(apply(lr_all,2,function(x){max(x,na.rm = TRUE)}),decreasing = TRUE)[nn]
  cutlevel2 <- max(apply(lr_all,1,function(x){sort(x,decreasing = TRUE)[nn]}))
  return(list(cut1=cutlevel1,cut2=cutlevel2))
}











