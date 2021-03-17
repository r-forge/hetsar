## from 'plm', /tool_misc.R

pres <- function(x) {  # pres.panelmodel
  ## extracts model residuals as pseries
  ## not necessary for plm models as residuals.plm returns a pseries,
  ## but used in residuals.pggls, residuals.pcce, residuals.pmg
  
  ## extract indices
  groupind <- attr(x$model, "index")[ , 1L]
  timeind  <- attr(x$model, "index")[ , 2L]
  
  # fix to allow operation with pggls, pmg
  # [TODO: one day, make this cleaner; with the describe framework?]
  if (!is.null(x$args$model))                 maybe_fd <- x$args$model
  if (!is.null(attr(x, "pmodel")$model.name)) maybe_fd <- attr(x, "pmodel")$model.name # this line is currently needed to detect pggls models
  
  ## Achim's fix: reduce id and time index to accommodate first-differences model's number of observations
  if(exists("maybe_fd") && maybe_fd == "fd") {
    groupi <- as.numeric(groupind)
    ## make vector =1 on first obs in each group, 0 elsewhere
    selector <- groupi - c(0, groupi[-length(groupi)])
    selector[1L] <- 1 # the first must always be 1
    ## eliminate first obs in time for each group
    groupind <- groupind[!selector]
    timeind <- timeind[!selector]
  }

  resdata <- data.frame(ee = x$residuals, ind = groupind, tind = timeind)
  pee <- pdata.frame(resdata, index = c("ind", "tind"))
  pres <- pee$ee
  return(pres)
}
