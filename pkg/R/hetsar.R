## wrapper for setting data and then calling fn_hsar()

hetsar <- function(formula, data, w, na.action,
                   index = NULL, ...) {
    #trace <- as.numeric(!quiet)
    cl <- match.call()
    
    ## input controls
    if (!is.matrix(w)) {
        if ("listw" %in% class(w)) {
            w <- listw2mat(w)
        }
        else {
            stop("w has to be either a 'matrix' or a 'listw' object")
        }
    }

    ## make clean data based on 'plm' interface
    ## (supports transformation functions in formula)

    ## estimate a pooling model by OLS (inexpensive)
    pmod <- plm(formula, data, index = index, model = "pooling")
    ## fetch data
    X <- model.matrix(pmod)
    y <- pmodel.response(pmod)
    ## fetch indices
    ind <- attr(pmod$model, "index")[, 1]
    tind <- attr(pmod$model, "index")[, 2]
    ## put data in "spatial" order as stack of cross-sections
    ## (ind is the "fast" index)
    oo <- order(tind, ind)
    X <- X[oo, , drop = FALSE]
    y <- y[oo]
    ind <- ind[oo]
    tind <- tind[oo]
    ## calc dimensions (might as well use pdim{plm})
    n <- length(unique(ind))
    k <- dim(X)[[2]]
    t <- max(tapply(X[, 1], ind, length))
    nT <- length(ind)

    ## check conformability with W
    if (dim(w)[[1]] != n) 
        stop("Non conformable spatial weights")

    ## check balanced
    balanced <- pdim(pmod)$balanced
    if (!balanced) { 
        stop("Estimation method unavailable for unbalanced panels")
        ## or whatever... for now let us assume balanced.
    }

    ## make data in fn_hsar format:

    ## make dataframe in correct order
    orddata <- data.frame(ind=ind, tind=tind, y=y)
    orddata <- cbind(orddata, X)
    
    ## define internal function
    panel2mat <- function(x) {
        tempmelt <- melt(cbind(orddata[,1:2], x),
                         id=c("ind","tind"))
        res <- cast(tempmelt, ind~tind)
        mres <- as.matrix(res)
        return(mres)
    }

    ## put y, X in matrix / list of matrices form
    m_y <- panel2mat(orddata$y)
    l_x <- lapply(orddata[, (4+1):(dim(orddata)[[2]]),
                          drop=FALSE],
                  FUN = panel2mat)

    ## FUTURE MOD:
    ## we assume intercept is in formula, if any,
    ## and we eliminate part in fn_hsar.r where it was
    ## added; now we take it away, see below;
    ## else the above line would
    ## be '4:(dim(orddata)[[2]]'

    ## make l_data for feeding to fn_hsar
    l_data <- list(m_y=m_y, l_x=l_x, m_W=w)   

    
#### this must work now:
#### fn_hsar(l_data)

    ## it does if the formula is w/o intercept!! because
    ## the intercept is added in fn_hsar.
    ## Check whether one could manage the intercept just
    ## like any other X, it would improve the consistency
    ## of code
    

    ## (init values module here)
    
    ## call estimator fun
    RES <- fn_hsar(l_data, ...) #, method, l_bounds, ...)

    ## make results object (take useful bits from spreml)

    #class(hetmod) <- "hetsar"
    #return(hetmod)

    ## for now, return the individual coefs' matrix:
    indcoef <- RES$theta

    names.beta <- names(coef(pmod))
    dimnames(indcoef)[[2]] <- c("SAR", names.beta[-1],
                                "Intercept", "variance")

    ## ...and warning if nonconvergence!

    ## return individual coefficents' variances too
    indcoefvar.std <- RES$standard
    indcoefvar.snd <- RES$sandwich
    
    ## simply return coefs:
    #return(indcoef)
    
###### results module, from pmg() ######
    ## make tcoef w/o variance parameter?
    tcoef <- t(indcoef)
    
    ## MG coefficients, omit variance (which comes last)
    coef <- rowMeans(tcoef)

    ## make MG variance (check if appropriate)
    coefmat <- array(dim = c(k+2, k+2, n))
    demcoef <- tcoef - coef
    for (i in 1:n) coefmat[, , i] <- outer(demcoef[, i],
                                           demcoef[, i])
    vcov <- apply(coefmat, 1:2, sum)/(n * (n - 1))
    
    residuals <- NULL #unlist(tres)
    df.residual <- NULL #nrow(X) - ncol(X)
    fitted.values <- NULL #y - residuals
    r2 <- NULL #1 - var(residuals)/var(y)*(T.-1)/(T.-k-1)
    #names(coef) <- rownames(vcov) <- colnames(vcov) <- coef.names
    #dimnames(tcoef) <- list(coef.names, id.names)
    
    #pmodel <- attr(pmod, "pmodel")
    #pdim <- attr(pmod, "pdim")
    #pmodel$model.name <- model.name
    hsarmod <- list(coefficients = coef,
                    residuals = residuals, 
                    fitted.values = fitted.values, vcov = vcov,
                    df.residual = df.residual, r.squared = r2,
                    model = model.frame(pmod),
                    sigma = NULL, indcoef = tcoef,
                    indcoefvar.std = indcoefvar.std,
                    indcoefvar.snd = indcoefvar.snd,
                    call = cl, pdim = pdim(pmod))
    #hsarmod <- structure(hsarmod, pmodel = pmodel)
    class(hsarmod) <- c("hetsar", "panelmodel")
    hsarmod
    
    return(hsarmod)
    
}



############################################

summary.hetsar <- function(object, ...){
  pmodel <- attr(object, "pmodel")
  #pdim <- attr(object, "pdim")
  std.err <- sqrt(diag(object$vcov))
  b <- object$coefficients
  z <- b/std.err
  p <- 2*pnorm(abs(z), lower.tail = FALSE)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable <- CoefTable
  y <- object$model[[1L]]
  object$tss <- sum(y^2) #tss(y)
  object$ssr <- NULL #sum(residuals(object)^2)
  object$rsqr <- NULL #1-object$ssr/object$tss
  class(object) <- c("summary.hetsar")
  return(object)
}

print.summary.hetsar <- function(x,
              digits = max(3, getOption("digits") - 2),
              width = getOption("width"), ...){
  pmodel <- attr(x, "pmodel")
  #pdim <- attr(x, "pdim")
#  formula <- pmodel$formula
  #model.name <- pmodel$model.name
  #cat(paste(model.hetsar.list[model.name], "\n", sep=""))
  cat("HSAR model")
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(x$pdim)
  cat("\nResiduals:\n")
  print(summary(unlist(residuals(x))))
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits)
  #cat(paste("Total Sum of Squares: ",
  #          signif(x$tss, digits), "\n", sep=""))
  #cat(paste("Residual Sum of Squares: ",
  #          signif(x$ssr, digits), "\n", sep=""))
  #cat(paste("Multiple R-squared: ",
  #          signif(x$rsqr, digits),"\n", sep=""))
  invisible(x)
}


residuals.hetsar <- function(object, ...) {
  return(pres(object))
}

