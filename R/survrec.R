"Survr" <-
function (id,time,event) 
{
 
    if(length(unique(id))!=length(event[event==0]))
      {
        stop("data doesn't match")
      }

    if(length(unique(event))>2 | max(event)!=1 | min(event)!=0)
      {
        stop("event must be 0-1")
      }

    ans<-cbind(id,time,event)

    class(ans) <- "Survr"
    invisible(ans)

}

"is.Survr" <-
function(x)
inherits(x, "Survr")


"psh.fit" <-
function(x,tvals)
{

      if (!is.Survr(x)) 
      {
           stop("\n x must be a Survr object")
      }

       n <- length(unique(x[,1]))
       failed<-c(x[,2][x[,3]==1])
       censored<-c(x[,2][x[,3]==0])
       m <- table(x[, 1])-1

       sfailed<-sort(failed)
       nfailed<-length(failed)

       summ <- .Fortran("distinctfailed",
                         as.integer(n),
                         as.integer(m),
                         as.double(failed),
                         as.double(sfailed),
                         as.integer(nfailed),
                         as.double(censored),
                         as.integer(0),
                         as.double(rep(0,nfailed)),
                         as.integer(rep(0,nfailed)),
                         as.integer(rep(0,n*nfailed)))
       numdistinct <- summ[[7]]
       distinct <- summ[[8]][1:numdistinct]
       numdeaths <- summ[[9]][1:numdistinct]
       vAtRisk <- summ[[10]][1:(n*numdistinct)]
       AtRisk<-matrix(vAtRisk,nrow=n,ncol=numdistinct)

       survfuncPSHple <- vector("numeric", numdistinct)
       AtRiskTotals <- t(AtRisk) %*% c(rep(1, n))

       survfuncPSHple.i<- 1-(numdeaths/AtRiskTotals)
       survfuncPSHple<-cumprod(survfuncPSHple.i)
        
       if(!missing(tvals))
        {
          tvalslen <- length(tvals)
          tvals.o<-sort(tvals)
          PSHpleAttvals<-surv.search(tvals.o,distinct,survfuncPSHple)
        }
       else
        {
          tvals<-NA
          PSHpleAttvals<-NA
        }

       ans<-list(n = n, m = m, failed = failed, censored = censored,
                 time = distinct, n.event=numdeaths, AtRisk = AtRisk,
                survfunc = survfuncPSHple, tvals = tvals, PSHpleAttvals
                 = PSHpleAttvals)
       class(ans)<-"survfitr"
       ans

}


"wc.fit" <-
function(x,tvals)
{
      if (!is.Survr(x)) 
      {
           stop("\n x must be a Survr object")
      }

       n <- length(unique(x[,1]))
       failed<-c(x[,2][x[,3]==1])
       censored<-c(x[,2][x[,3]==0])
       m <- table(x[, 1])-1
   
       sfailed<-sort(failed)
       nfailed<-length(failed)

        summ <- .Fortran("distinctfailed",
                         as.integer(n),
                         as.integer(m),
                         as.double(failed),
                         as.double(sfailed),
                         as.integer(nfailed),
                         as.double(censored),
                         as.integer(0),
                         as.double(rep(0,nfailed)),
                         as.integer(rep(0,nfailed)),
                         as.integer(rep(0,n*nfailed)))

        numdistinct <- summ[[7]]
        distinct <- summ[[8]][1:numdistinct]
        numdeaths <- summ[[9]][1:numdistinct]
        vAtRisk <- summ[[10]][1:(n*numdistinct)]
        AtRisk<-matrix(vAtRisk,nrow=n,ncol=numdistinct)

        wcPLE<- .Fortran("wcple",
                       as.integer(n),
                       as.integer(m),
                       as.double(failed),
                       as.integer(nfailed),
                       as.double(censored),
                       as.integer(numdistinct),
                       as.double(distinct),
                       as.integer(c(vAtRisk)),
                       as.double(rep(0,n*numdistinct)),
                       as.double(rep(0,n*numdistinct)),
                       as.double(rep(0,n)))

        dstar<-matrix(wcPLE[[9]],n,numdistinct)
        rstar<-matrix(wcPLE[[10]],n,numdistinct)
        mstar<-wcPLE[[11]]
       
        dstarcum <- t(dstar) %*% (c(rep(1, n))/mstar)
        rstarcum <- t(rstar) %*% (c(rep(1, n))/mstar)

        survfuncWCple.i<-function(i,x,y)
                          {
                            exp(sum(log(1 - (x[1:i]/y[1:i]))))
                          }

        survfuncWCple<-sapply(1:numdistinct,survfuncWCple.i,x=dstarcum,y=rstarcum)
             
        if(!missing(tvals))
         { 
          tvalslen <- length(tvals)
          tvals.o<-sort(tvals)
          WCpleAttvals<-surv.search(tvals.o,distinct,survfuncWCple)
         }
        else
         {
          tvals<-NA
          WCpleAttvals<-NA     
         } 

        ans<-list(n = n, m = m, failed = failed, censored = censored, 
                 time = distinct, n.event=numdeaths, AtRisk = AtRisk,
               survfunc = survfuncWCple, tvals = tvals, WCpleAttvals = 
                   WCpleAttvals)
        class(ans)<-"survfitr"
        ans
}


"mlefrailty.fit"<-
function (x, tvals, lambda = NULL, alpha = NULL, alpha.min, alpha.max, 
    tol = 1e-07, maxiter = 500, alpha.console=TRUE) 
{
    if (!is.Survr(x)) {
        stop("\n x must be a Survr object")
    }
    n <- length(unique(x[, 1]))
    failed <- c(x[, 2][x[, 3] == 1])
    censored <- c(x[, 2][x[, 3] == 0])
    m <- table(x[, 1])-1
    sfailed <- sort(failed)
    nfailed <- length(failed)
    summ <- .Fortran("distinctfailed", as.integer(n), as.integer(m), 
        as.double(failed), as.double(sfailed), as.integer(nfailed), 
        as.double(censored), as.integer(0), as.double(rep(0, 
            nfailed)), as.integer(rep(0, nfailed)), as.integer(rep(0, 
            n * nfailed)))
    numdistinct <- summ[[7]]
    distinct <- summ[[8]][1:numdistinct]
    numdeaths <- summ[[9]][1:numdistinct]
    vAtRisk <- summ[[10]][1:(n * numdistinct)]
    AtRisk <- matrix(vAtRisk, nrow = n, ncol = numdistinct)
    if (is.null(lambda[1])) {
        lambda <- as.vector(numdeaths/apply(AtRisk,2,sum))
    }
    if (is.null(alpha)) {
 
        if(alpha.console)
              cat("\nNeeds to Determine a Seed Value for Alpha")
        if (missing(alpha.min)) {
            alpha.min <- 0.5
        }
        if (missing(alpha.max)) {
            alpha.max <- max(distinct)
        }
        tol.max <- (alpha.max - alpha.min)/50
        Seed <- .Fortran("searchforseed", as.integer(n), as.integer(m), 
            as.integer(numdistinct), as.double(distinct), as.integer(numdeaths), 
            as.integer(AtRisk), as.double(lambda), as.double(alpha.min), 
            as.double(alpha.max), as.double(tol.max), as.double(0), 
            as.integer(0))
        alpha <- Seed[[11]]
        if(alpha.console)
               cat("\n Seed Alpha: ", alpha)
        IER <- Seed[[12]]
        if (IER != 0 && alpha.console) {
            warning("Problem with the seed value for alpha")
            if (IER == 129) {
                warning("bL is greater than or equal to bR. Minimum as bL")
            }
            if (IER == 130) {
                warning("tol is greater than the interval bL to bR")
            }
            if (IER == 131) {
                warning("the function is not unimodal. Check your results")
            }
        }
    }
    alphadel <- alpha/4
    alphaseeds <- c(alpha, alpha - alphadel, alpha - 2 * alphadel, 
        alpha - 3 * alphadel, alpha + alphadel, alpha + 2 * alphadel, 
        alpha + 3 * alphadel)
    status <- 0
    ind <- 0
    while ((status == 0) && (ind < 7)) {
        ind <- ind + 1
        alpha <- alphaseeds[ind]
        Estimates <- .Fortran("emalgo", as.integer(n), as.integer(m), 
            as.integer(numdistinct), as.double(distinct), as.integer(numdeaths), 
            as.integer(AtRisk), as.double(lambda), as.double(alpha), 
            as.double(tol), as.integer(maxiter), as.integer(status))
        status <- Estimates[[11]]
    }
    alpha <- Estimates[[8]]
    if (alpha.console)
      { 
        cat("\n ")
        cat("\n Alpha estimate=", alpha)
        cat("\n ")
      }
    lambda <- Estimates[[7]]
    if (!missing(tvals)) 
        tvalslen <- length(tvals)
    if (!(status == 1)) {
        cat("\n\n WARNING: No estimates will be provided!")
        cat("\n Value of (status,alpha) from iteration is ", 
            c(status, alpha), "\n\n")
        alpha <- NA
        survfuncMLE <- c(rep(NA, numdistinct))
        if (!missing(tvals)) 
            MLEAttvals <- c(rep(NA, tvalslen))
        else {
            tvals <- NA
            MLEAttvals <- NA
        }
    }
    else {
        if (alpha >= 1e+05) 
            alpha <- 1e+05
        temp <- .Fortran("mlevalue", as.integer(numdistinct), 
            as.double(alpha), as.double(lambda), as.double(rep(0, 
                numdistinct)))
        survfuncMLE <- temp[[4]]
        if (!missing(tvals)) {
            tvalslen <- length(tvals)
            tvals.o <- sort(tvals)
            MLEAttvals <- surv.search(tvals.o, distinct, survfuncMLE)
        }
        else {
            tvals <- NA
            MLEAttvals <- NA
        }
    }
    ans <- list(n = n, m = m, failed = failed, censored = censored, 
        time = distinct, n.event = numdeaths, AtRisk = AtRisk, 
        status = status, alpha = alpha, lambda=lambda, survfunc = survfuncMLE, 
        tvals = tvals, MLEAttvals = MLEAttvals)
    class(ans) <- "survfitr"
    ans
}




"surv.search" <-
function (tvals,time,surv)
{
  time.c<-c(0,time)
  surv.c<-c(1,surv)
  mm<-outer(tvals,time.c,">=")
  pos<-apply(mm,1,sum)
  return(surv.c[pos])
}


"survfitr" <-
function (formula, data, type="MLEfrailty",...) 
{
   method <- charmatch(type, c("pena-strawderman-hollander", 
             "wang-chang", "MLEfrailty"), nomatch= 0)
   if(method == 0)
	{
	  stop("estimator must be pena-strawderman-hollander wang-chang or mlefrailty")
        }

    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) || 
        inherits(formula, "Survr")) {

    stop("formula.default(object): invalid formula")
     }

    m <- match.call(expand = FALSE)
    m$type<- m$... <- NULL
    Terms <- terms(formula, "strata")
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)
    Y <- model.extract(m, response)
    if (!is.Survr(Y)) 
        stop("Response must be a survival recurrent object")
    ll <- attr(Terms, "term.labels")
    group <- m[ll][, 1]


    if(method==1) FUN<-psh.fit
    if(method==2) FUN<-wc.fit
    if(method==3) FUN<-mlefrailty.fit

    if (!is.null(group)) {
        k <- levels(group)
        ans <- NULL
        for (i in 1:length(k)) {
            temp <- Y[group == k[i], ]
            temp1 <- Survr(temp[, 1], temp[, 2], temp[, 3])
            ans[[i]] <- FUN(temp1,...)
        }
        names(ans) <- k
        class(ans) <- "survfitr"
        attr(ans, "strata") <- length(k)
        attr(ans, "group") <- ll
    }
    else {
        temp<-Survr(Y[,1],Y[,2],Y[,3]) 
        ans <- FUN(temp,...)
    }
    ans
}





"plot.survfitr" <-
function (x, prob=FALSE,...) 
{
    dostep <- function(x, y) {
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }

y.lab<-ifelse(prob,"Probability Estimates","Survivor Probability Estimates")


if(!prob){
    if (!is.null(attr(x, "strata")))  {
        plot(dostep(x[[1]]$time, x[[1]]$surv), type = "n", xlab = "Time", 
            ylab = y.lab, ...)
        for (i in 1:attr(x, "strata")) {
            lines(dostep(x[[i]]$time, x[[i]]$surv),lty=i)
        }
    }
    else {
        y <- x$survfunc
        plot(dostep(x$time, y), type = "l", ylim = c(0, max(y)), 
            xlab = "Time", ylab =y.lab)
    }

}

else{
    if (!is.null(attr(x, "strata")))  {
        plot(dostep(x[[1]]$time, 1-x[[1]]$surv), type = "n", xlab = "Time", 
            ylab = y.lab, ...)
        for (i in 1:attr(x, "strata")) {
            lines(dostep(x[[i]]$time, 1-x[[i]]$surv),lty=i)
        }
    }
    else {
        y <- x$survfunc
        plot(dostep(x$time, 1-y), type = "l", ylim = c(0, max(y)), 
            xlab = "Time", ylab = y.lab)
    }

}

    return(invisible())
}


"lines.survfitr"<-
function (x, prob=FALSE, ...) 
{
    dostep <- function(x, y) {
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }
if(!prob)
    lines(dostep(x$time, x$survfunc), ...)
else
    lines(dostep(x$time, 1-x$survfunc), ...)
    return(invisible())
}



"print.survfitr" <-
function (x,scale=1,digits = max(options()$digits - 4, 3), ...) 
{

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  plab<-c("n","events","mean","se(mean)","median","recurrences: min","max","median")

  pfun <- function(x)
    #compute the mean, se(mean) and median survival
	{
          minmin <- function(y, xx)
            {
	     if(any(!is.na(y) & y == 0.5)) {
		if(any(!is.na(y) & y < 0.5))
		  0.5 * (min(xx[!is.na(y) & y == 0.5]) + min(xx[!is.na(y) & y < 0.5]))
		else 0.5 * (min(xx[!is.na(y) & y == 0.5]) + max(xx[!is.na(y) & y == 0.5]))
			}
			else min(xx[!is.na(y) & y <= 0.5])
		}

                stime<-x$time/scale
           	    n <- length(stime)
          	    n.risk<-apply(x$AtRisk,2,sum)
                hh <- c(x$n.event[ - n]/(n.risk[ - n] * (n.risk[ - n] - x$n.event[ - n])), 0)
	       
                med <- minmin(x$survfunc, x$time)
               

                dif.time <- c(diff(c(0, stime)), 0)
                mean <- dif.time * c(1, x$survfunc)

                varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)                 

                ans<-c(x$n,sum(x$m),sum(mean),sqrt(varmean),med,min(x$m),max(x$m),median(x$m))
                ans
}


if(is.null(attr(x,"strata")))
  {
   # no strata 
    x1<-rbind(pfun(x))
    cat("Survival for recurrent event data")
    cat("\n")
    dimnames(x1)<-list(" ",plab)
    print(x1)
    cat("\n")
  }

else
  {
      cat("Survival for recurrent event data. Group=" ,attr(x,"group"))
      cat("\n")  
      x1<-NULL
      for (i in 1:attr(x,"strata"))
        {
         temp<-pfun(x[[i]])
         x1<-rbind(x1,temp)
        }
      temp<-names(x)
      dimnames(x1)<-list(temp,plab)
      print(x1)
      cat("\n")
   }

invisible(x)

}


"summary.survfitr" <-
function (object,...) 
{
  x<-object
  if (!inherits(x, "survfitr")) 
        stop("Invalid data")
  
  plab<-c("time","n.event","n.risk","surv")

  if(!is.null(attr(x,"strata")))
     {
       ans<-list(NA)
       for (i in 1:attr(x,"strata"))
         {
          n.risk<-apply(x[[i]]$AtRisk,2,sum)
          temp<-cbind(x[[i]]$time,x[[i]]$n.event,n.risk,x[[i]]$surv)
          dimnames(temp)<-list(rep("",nrow(temp)),plab)
          ans[[i]]<-temp
         }
       names(ans)<-names(x)
       class(ans)<-"summary.survfitr"
       attr(ans,"strata")<-attr(x,"strata")
     }
  else
     {
          n.risk<-apply(x$AtRisk,2,sum)
          temp<-cbind(x$time,x$n.event,n.risk,x$surv)
          dimnames(temp)<-list(rep("",nrow(temp)),plab)
          ans<-temp
     }  
 ans
}



"print.summary.survfitr" <-
function (x,scale=1,digits = max(options()$digits - 4, 3), ...) 
{
  savedig <- options(digits = digits)
  on.exit(options(savedig))

if(is.null(attr(x,"strata")))
  {
   # no strata 
    cat("\n")
     temp<-x
    class(temp)<-NULL
    print(temp)
    cat("\n")
  }

else
 { 
  for (i in 1:attr(x,"strata"))
    {
     cat("\n      Group=",names(x)[i])
     cat("\n")
     print(x[[i]])
     cat("\n")
    }
 }
 invisible(x)
}


############ First.lib ###############

.First.lib <- function(lib, pkg){
   library.dynam("survrec", pkg, lib)
   cat("   Survival analysis for recurrent event data installed\n")
   cat("   created by Juan R Gonzalez, Edsel A Peña and Robert L Strawderman\n")
}
############ End of .First.lib ###############
