#' @export
PGEE<-PGEE_own<- function(formula, id, data, na.action = NULL, family = gaussian(link = "identity"),
           corstr = "independence", Mv = NULL, beta_int = NULL, R = NULL, scale.fix = TRUE,
           scale.value = 1, lambda, pindex = NULL, eps = 10^-6, maxiter = 30, tol = 10^-3, 
           silent = TRUE)  {
  #library(tictoc)
    #pass the call
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    
    #initilize non-arguments
    m$beta_int <- m$family <- m$link <- m$varfun<-
      m$corstr <- m$Mv<- m$R <-
      m$scale.fix <- m$scale.value <-
      m$lambda <-m$eps<-m$pindex<-
      m$maxiter <- m$tol <-m$silent <- NULL
    
    #initilize cluster
    if(is.null(m$id)) m$id<-as.name("id")
    
    #throw out warning on na.omit
    if(!is.null(m$na.action) && m$na.action != "na.omit") {
      warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
      m$na.action <- as.name("na.omit")
    }
    
    
    #assign arguments
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    #extract data
    y <- model.extract(m, "response")
    X<- model.matrix(Terms, m)
    #extract cluster
    id<-model.extract(m, id)
    
    #cluster id is required and must be in same length
    if(is.null(id)) {
      stop("Id variable not found!")
    }
    if(length(id) != length(y))  stop("Id and y do not have the same length!")
    
    if(!(is.double(X)))  X <- as.double(X)
    if(!(is.double(y)))  y <- as.double(y)
    if(!(is.double(id))) id <- as.double(id)
    
    #K as the variable number adjusted by the presence of intercept
    N<-length(unique(id))
    if (colnames(X)[1]== "(Intercept)") K=dim(X)[2]-1 else K=dim(X)[2]
    nx=ncol(X)
    
    #avec as size of each cluster
    avec <- as.integer(unlist(lapply(split(id, id), "length")))
    #largest cluster
    maxclsz <-max(avec)
    maxcl <- maxclsz
    nt<-avec
    #total obs number
    nobs<-sum(nt)
    
    #reset x names
    xnames <- dimnames(X)[[2]]
    if(is.null(xnames) && colnames(X)[1]=="(Intercept)") {
      xnames <- paste("x", 0:K, sep = "")
      dimnames(X) <- list(NULL, xnames)
    } else
      if(is.null(xnames) && colnames(X)[1]!="(Intercept)") {
        xnames <- paste("x", 1:K, sep = "")
        dimnames(X) <- list(NULL, xnames)
      }
    
    
    if(!(is.double(N)))      N <- as.double(N)
    if(!(is.double(maxcl)))  maxcl <- as.double(maxcl)
    if(!(is.double(nobs)))   nobs <- as.double(nobs)
    
    #set default if missing; validate each arguments
    if(missing(lambda)) stop("A value is not assiged for lambda!")
    
    if(missing(pindex)) pindex=NULL
    
    if(missing(family)) family=gaussian(link="identity")
    
    if(missing(corstr)) corstr="independence"
    
    if(missing(Mv)) Mv<-NULL
    
    if(corstr=="stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'stat_M_dep' but Mv is not specified!")
    
    if(corstr=="non_stat_M_dep" && is.null(Mv)) stop("corstr is assumed to be 'non_stat_M_dep' but Mv is not specified!")
    
    if((corstr!="stat_M_dep" && corstr!="non_stat_M_dep") && !is.null(Mv))  stop("Mv is specified while corstr is assumed to be neither 
'stat_M_dep' nor 'non_stat_M_dep'!")
    
    if(corstr=="non_stat_M_dep" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'non_stat_M_dep' for unbalanced data!")
    
    if(corstr=="unstructured" && length(unique(nt)) !=1) stop("corstr cannot be assumed to be 'unstructured' for unbalanced data!")
    
    if(missing(R)) R<-NULL
    
    if(corstr=="fixed" && is.null(R))  stop("corstr is assumed to be 'fixed' but R is not specified!")
    if(corstr!="fixed" && !is.null(R)) stop("R is specified although corstr is not assumed to be 'fixed'!")
    
    if(!is.null(R)) {
      Rr <- nrow(R)
      if(Rr != ncol(R)) stop("R is not square!")
      if(Rr < maxclsz)  {stop("R is not big enough to accommodate some clusters!")} else
        if(Rr > maxclsz)  {stop("R is larger than the maximum cluster!")}
    }
    
    if(missing(scale.fix))  scale.fix <- TRUE
    scale.fix <- as.integer(scale.fix)
    
    if(missing(scale.value)) scale.value=1
    scale.value<-as.integer(scale.value)
    
    if(missing(eps)) eps=10^-6
    eps<-as.double(eps)
    
    if(missing(maxiter)) maxiter<-30
    maxiter<-as.integer(maxiter)
    
    if(missing(tol))  tol=10^-3
    tol=as.double(tol)
    
    if(missing(silent))  silent <-TRUE
    silent<-as.integer(silent)
    
    if (is.character(family)) family <- get(family)
    if (is.function(family))  family <- family()
    
    links <- c("identity","log","logit","inverse","probit","cloglog")
    fams <- c("gaussian","poisson","binomial","Gamma","quasi")
    varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
    corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", 
                 "AR-1", "unstructured")
    
    linkv <- as.integer(match(c(family$link), links, -1))
    if(linkv < 1) stop("unknown link!")
    
    famv <- match(family$family, fams, -1)
    if(famv < 1) stop("unknown family")
    if(famv <= 4) varfunv <- famv
    else varfunv <- match(family$varfun, varfuns, -1)
    if(varfunv < 1) stop("unknown varfun!")
    
    corstrv <- as.integer(match(corstr, corstrs, -1))
    if(corstrv < 1) stop("unknown corstr!")
    
    Mv <- as.integer(Mv)
    
    if (!is.null(beta_int))
    {
      beta <- matrix(beta_int, ncol = 1)
      if(nrow(beta) != nx) {stop("Dimension of beta != ncol(X)!")}
      #message("user\'s initial regression estimate")
      
    }
    else {
      #message("running glm to get initial regression estimate!")
      ### <tsl>	beta <- as.numeric(glm(m, family = family)$coef)
      mm <- match.call(expand.dots = FALSE)
      mm$R <- mm$beta_int <- mm$tol <- mm$maxiter <- mm$link <- 
        mm$varfun <-mm$corstr <- mm$Mv <- mm$silent <-mm$scale.fix <- 
        mm$scale.value <- mm$id<-
        mm$lambda <-mm$pindex<-mm$eps<-NULL
      mm[[1]]<-as.name("glm")
      beta <- eval(mm, parent.frame())$coef
      ### </tsl>
      print(beta)
      
    }
    
    beta_int=matrix(beta, ncol = 1)
    beta_new<-beta_int
    

    # initial estimate of working correlation matrix
    R.fi.hat=mycor_gee2(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R=R,scale.fix=scale.fix,scale.value=scale.fix)
    Rhat=R.fi.hat$Ehat
    fihat=R.fi.hat$fi
    
    
    
    # initial estimates PGEE
    S.H.E.val=S_H_E_M(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps)
    S<-S.H.E.val$S
    H<-S.H.E.val$H
    E<-S.H.E.val$E
    
    diff<-1
    iter<-0
    # interation to update beta, Newton-Raphson algorithm

    
    while(iter < maxiter) {
      beta_old<-beta_new
      #update beta by equation 5.1
      #tic("update beta")
      beta_new<-matrix(beta_old)+geninv(H+N*E)%*%(S-N*E%*%matrix(beta_old))
      #estimate R with new beta p1
      ##toc()
      #tic("working matrix")
      R.fi.hat=mycor_gee2(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R,scale.fix,scale.value)
      Rhat=R.fi.hat$Ehat
      fihat=R.fi.hat$fi
      #toc()
      #estimate SHEM with new beta p2
      
      S.H.E.M.val=S_H_E_M(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps)
      S<-S.H.E.M.val$S
      H<-S.H.E.M.val$H
      E<-S.H.E.M.val$E
      M<-S.H.E.M.val$M
      
      #set difference
      diff<-sum(abs(beta_old-beta_new)) 
      
      iter<-iter+1
      if (silent==0) cat("iter",iter,"beta_new",beta_new,"diff",diff,"\n")
      #break if difference is small enough
      if (diff <= tol) break
    } #end of while
    
    #toc()
    

    estb=beta_new
    nv=naive.var<-geninv(H+N*E)
    rv=robust.var<-geninv(H+N*E)%*%M%*%geninv(H+N*E)
    final_iter=iter
    final_diff=diff
    
    fit <- list()
    attr(fit, "class") <- c("PGEE","gee","glm")
    fit$title <- "PGEE: PENALIZED GENERALIZED ESTIMATING EQUATIONS FOR LONGITUDINAL DATA"
    fit$version <- "Version: 1.5"
    links <- c("Identity", "Logarithm", "Logit", "Reciprocal", "Probit","Cloglog")
    varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
    corstrs <- c("Independent", "Fixed", "Stationary M-dependent",
                 "Non-Stationary M-dependent", "Exchangeable", "AR-1",
                 "Unstructured")
    fit$model <- list()
    fit$model$link <- links[linkv]
    fit$model$varfun <- varfuns[varfunv]
    fit$model$corstr <- corstrs[corstrv]
    if(!is.na(match(c(corstrv), c(3, 4))))
      fit$model$M <- Mv
    fit$call <- call
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    #fit$contrasts <- attr(X, "contrasts")
    fit$nobs <- nobs
    fit$iterations <- final_iter
    fit$coefficients <- as.vector(estb)
    fit$nas <- is.na(fit$coefficients)
    names(fit$coefficients) <- xnames
    eta <- as.vector(X %*% fit$coefficients)
    fit$linear.predictors <- eta
    ##Rchange
    mu <- as.vector(family$linkinv(eta))
    ##
    fit$fitted.values <- mu
    fit$residuals <- y - mu
    fit$family <- family
    fit$y <- as.vector(y)
    fit$id <- as.vector(id)
    fit$max.id <- maxcl
    fit$working.correlation <- Rhat[1:maxclsz,1:maxclsz,which(avec==maxclsz)[1]]
    fit$scale <- fihat
    fit$epsilon<-eps
    fit$lambda.value<-lambda
    fit$robust.variance <- rv
    fit$naive.variance <- nv
    fit$xnames <- xnames
    fit$error <- final_diff
    dimnames(fit$robust.variance) <- list(xnames, xnames)
    dimnames(fit$naive.variance) <- list(xnames, xnames)
    #toc()
    fit
    
}

#working correlation matrix against the structure
mycor_gee2 <-function(N,nt,y,X,family,beta_new,corstr,Mv,maxclsz,R,scale.fix,scale.value) {
  eta=X%*%beta_new
  mu=family$linkinv(eta)
  sd=sqrt(family$variance(mu))
  res<-(as.vector(y)-mu)/sd 
  
  if(scale.fix==0) {
    fi<-sum(res^2)/(sum(nt))
  } else
    if(scale.fix==1) {
      fi<-scale.value
    }
  
  aindex=cumsum(nt)
  index=c(0,aindex[-length(aindex)])
  
  if (corstr=="independence") 
  {alfa_hat<-0} else 
    if (corstr=="exchangeable") 
    {
      sum1<-0
      sum3<-0
      for ( i in  1:N)          {
        for ( j in  1:nt[i])      {    
          for ( jj in 1:nt[i])      {    
            if  ( j!=jj)              {
              #cat("i",i,"j",j,"jj",jj,"\n")
              sum2<-res[j+index[i]]*res[jj+index[i]]
              #cat("i",i,"j",j,"jj",jj,"\n")
              sum1<-sum1+sum2
              #cat("i",i,"j",j,"jj",jj,"sum2",sum2,"sum1",sum1,"\n")
            }
          }
        }
        sum4<-nt[i]*(nt[i]-1)
        sum3<-sum3+sum4
      } #i
      alfa_hat<-sum1/(sum3*fi)
    } else
      if (corstr=="AR-1") 
      { 
        sum5<-0
        sum6<-0
        for ( i in  1:N)           {
          for ( j in  1:nt[i])       {  
            for ( jj in 1:nt[i])       {  
              if( j>jj && abs(j-jj)==1)  {
                #cat("i",i,"j",j,"jj",jj,"\n")
                sum7<-res[j+index[i]]*res[jj+index[i]]
                sum5<-sum5+sum7           
                #cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
              }
            }
          }
          sum8<-(nt[i]-1)
          sum6<-sum6+sum8
        } #i
        alfa_hat<-sum5/(sum6*fi)
      } else
        if (corstr=="stat_M_dep") 
        {  
          alfa_hat=matrix(0,Mv,1)
          for(m in 1:Mv) {
            sum12<-0
            sum14<-0
            for ( i in  1:N)           {
              for ( j in  1:nt[i])       {  
                for ( jj in 1:nt[i])       {  
                  if( j>jj && abs(j-jj)==m)  {
                    #cat("m",m,"i",i,","j",j,"jj",jj,"\n") 
                    sum11<-res[j+index[i]]*res[jj+index[i]]
                    sum12<-sum12+sum11 
                    #cat("m",m,"i",i,"j",j,"jj",jj,"sum11",sum11,"sum12", sum12, "\n")          
                  } #if
                }
              }
              sum13<-nt[i]-1
              sum14<-sum14+sum13
            } #i
            alfa_hat[m]<-sum12/(sum14*fi) 
          } #m
        }  else
          if (corstr=="non_stat_M_dep") 
          {  
            alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
            for( m in 1:Mv)            {
              for ( j in  1:nt[1])       {  
                for ( jj in 1:nt[1])       {  
                  if( j>jj && abs(j-jj)==m)  { 
                    sum16<-0                 
                    for ( i  in 1:N)     {
                      #cat("m",m,"j",j,"jj",jj,"i",i"\n") 
                      sum15<-res[j+index[i]]*res[jj+index[i]]
                      sum16<-sum15+sum16          
                      #cat("j",j,"jj",jj,"i",i,"sum15",sum15,"sum16",sum16,"\n")
                    } #i
                    #cat("j",j,"jj",jj,"sum16",sum16,"\n")
                    alfa_hat[j,jj]<-sum16/(N*fi) 
                  }
                }
              }
            }
          } else
            if (corstr=="unstructured") 
            {  
              alfa_hat<-matrix(0,nt[1],nt[1]) #not allowed for unequal number of cluster sizes.
              for ( j in 1:nt[1])  {  
                for ( jj in 1:nt[1]) {
                  sum20<-0                
                  if (j > jj)          {
                    for ( i  in 1:N )    {
                      #cat("i",i,"j",j,"jj",jj,"\n") 
                      sum21<-res[j+index[i]]*res[jj+index[i]]
                      sum20<-sum21+sum20           
                    } #i
                    #cat("j",j,"jj",jj,"sum20",sum20,"\n")
                    alfa_hat[j,jj]<-sum20/(N*fi) 
                  }
                }
              }
            } else
              if (corstr=="fixed")
              {alfa_hat=NULL
              }
  
  Ehat<-array(0,c(maxclsz,maxclsz,N))
  
  for(i in 1:N){
    cor1<-matrix(0,nt[i],nt[i])
    if (corstr=="independence")                                        
    {cor1<-diag(nt[i])} else
      if (corstr=="exchangeable")                                        
      { for (t1 in 1:nt[i]) {
        for (t2 in 1:nt[i]) {
          if (t1!=t2) 
          {cor1[t1,t2]<-alfa_hat} else 
          {cor1[t1,t2]<-1}
        }
      }
      } else
        if (corstr=="AR-1")                                      
        { for (t1 in 1:nt[i]) {
          for (t2 in 1:nt[i]) {
            cor1[t1,t2]<-alfa_hat^abs(t1-t2)   
          }
        }
        }  else
          if (corstr=="stat_M_dep")                                     
          { for (t1 in 1:nt[i]) {
            for (t2 in 1:nt[i]) {
              if (abs(t1-t2)==0)
              {cor1[t1,t2]<-1} else
                for(m in 1:Mv) {
                  if (abs(t1-t2)==m)
                  {cor1[t1,t2]<-alfa_hat[m]} 
                }
            }
          }
          } else
            if (corstr=="non_stat_M_dep")                                     
            { 
              cor1=alfa_hat+t(alfa_hat)
              diag(cor1)=1
            } else
              if (corstr=="unstructured")                                        
              {cor1=alfa_hat+t(alfa_hat)
              diag(cor1)=1
              } else
                if (corstr=="fixed")
                {cor1=R
                }
    
    Ehat[1:nt[i],1:nt[i],i]<-cor1 
  }
  
  return(list("Ehat"=Ehat,"fi"=fi))
  
}

#solve for H and E, with SCAD.
# S is GEE with the estimated R ; E is penalty matrix;
S_H_E_M <- function(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps=10^-6) {
  #tic("whole SHEM")
  aindex=cumsum(nt)
  index=c(0,aindex[-length(aindex)])
  
  eta=X%*%beta_new
  mu=family$linkinv(eta)
  
  #This is E on Wang et al.(2012) eq 5.2
  E1<-diag(q_scad(abs(as.vector(beta_new)),lambda)/(abs(as.vector(beta_new))+eps))
  
  #if some variable are not to be penalized
  if(is.null(pindex)==TRUE) {
    E<-E1 
  } else 
    if(is.null(pindex)!=TRUE) {
      E1[,pindex]<-0
      E<-E1
    }
  
  sum201<-matrix(0,nx,1)      #gradient:S
  sum301<-matrix(0,nx,nx)     #naive variance:H
  sum401<-matrix(0,nx,nx)     #a component for robust variance:M
  #
  for (i in 1:N) {
    ym<-matrix(0,nt[i],1)
    bigD<-matrix(0,nt[i],nx)
    bigA<-matrix(0,nt[i],nt[i])
    for (j in 1:nt[i]) {
      #cat("j",j,"\n")
      ym[j]<- y[j+index[i]]-mu[j+index[i]] 
      bigA[j,j]<-family$variance(mu)[j+index[i]]
      for (k in 1:nx) {
        bigD[j,k]<-family$mu.eta(eta)[j+index[i]]*X[j+index[i],k]
        #cat("i",i,"j",j,"k",k,"\n")
      } # for k
    } # for j
    ##working covariance matrix
    bigV<-sqrt(bigA)%*%Rhat[1:nt[i],1:nt[i],i]%*%sqrt(bigA)
    #bigV<-fihat*bigV
    
    ##This is S in Wang et al.(2012) eq4
    sum200<-t(bigD)%*%geninv(bigV)%*%ym      #this is like gradient
    sum201<-sum201+sum200
    
    ##This is H in Wang et al.(2012)
    sum300<-t(bigD)%*%geninv(bigV)%*%bigD    #this is for information matrix.
    sum301<-sum301+sum300
    
    ##Speed up the code##
    SSA=sqrt(geninv(bigA))
    SRhat=geninv(Rhat[1:nt[i],1:nt[i],i])
    SSAym=(SSA%*%ym)
    
    sum400<-t(bigD)%*%SSA%*%SRhat%*%(SSAym%*%t(SSAym))%*%SRhat%*%SSA%*%bigD
    sum401<-sum401+sum400
    #cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
  }  #end of i

  S<-fihat*sum201
  H<-fihat*sum301
  E<-E
  M<-fihat*sum401
  #toc()
  return(list("S"=S,"H"=H,"E"=E,"M"=M))
}

#SCAD penalty
q_scad <-function(theta,lambda,a=3.7)
{
  #length of parameter
  p<-length(theta)
  #get the absolute value
  theta<-abs(theta)
  #create vector of zeros
  b1<-rep(0,p)
  #if theta is greater then lambda set it to 1
  b1[theta>lambda]<-1
  #create an another vector of zeros
  b2<-rep(0,p)
  #if theta is less than a*lambda, set it to 1.
  b2[theta<(lambda*a)]<-1
  lambda*(1-b1)+((lambda*a)-theta)*b2/(a-1)*b1
}

#' @export
summary.PGEE <- summary.PGee <- function(object, correlation = TRUE, ...)
{
  coef <- object$coefficients
  resid <- object$residuals
  n <- length(resid)
  p <- object$rank
  if(is.null(p))
    p <- sum(!is.na(coef))
  if(!p) {
    warning("This model has zero rank --- no summary is provided")
    return(object)
  }
  nas <- is.na(coef)
  cnames <- names(coef[!nas])
  coef <- matrix(rep(coef[!nas], 5), ncol = 5)
  dimnames(coef) <- list(cnames, c("Estimate",
                                   "Naive S.E.",  "Naive z",
                                   "Robust S.E.", "Robust z"))
  rse <- sqrt(diag(object$robust.variance))
  nse <- sqrt(diag(object$naive.variance))
  coef[,2] <- nse
  coef[,3] <- coef[,1]/coef[,2]
  coef[,4] <- rse
  coef[,5] <- coef[,1]/coef[,4]
  summary <- list()
  summary$call <- object$call
  summary$version <- object$version
  summary$nobs <- object$nobs
  summary$residual.summary <- quantile(as.vector(object$residuals))
  names(summary$residual.summary) <- c("Min", "1Q", "Median", "3Q", "Max")
  summary$model<- object$model
  summary$title <- object$title
  summary$coefficients <- coef
  summary$working.correlation <- object$working.correlation
  summary$scale <- object$scale
  summary$lambda.value <-object$lambda.value
  summary$error <- paste("Error code was", object$error)
  summary$working.correlation <- object$working.correlation
  summary$iterations <- object$iterations
  if ( correlation ) {
    ##	rob.var <- object$robust.variance
    ##	nai.var <- object$naive.variance
    ##	summary$robust.correlation <- rob.var /
    ##	outer(sqrt(diag(rob.var)),sqrt(diag(rob.var)))
    ##	dimnames(summary$robust.correlation) <- list(object$xnames,object$xnames)
    ##	summary$naive.correlation <- nai.var /
    ##	outer(sqrt(diag(nai.var)),sqrt(diag(nai.var)))
    ##	dimnames(summary$naive.correlation) <- list(object$xnames,object$xnames)
  }
  attr(summary,"class") <- "summary.PGEE"
  summary
}
