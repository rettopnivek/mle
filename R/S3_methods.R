#------------#
# S3 methods #
#------------#

# Index
# Lookup - 01:  print.mle
# Lookup - 02:  coef.mle
# Lookup - 03:  sd.mle
# Lookup - 04:  vcov.mle
# Lookup - 05:  confint.mle
# Lookup - 06:  logLik.mle
# Lookup - 07:  deviance.mle
# Lookup - 08:  AIC.mle
# Lookup - 09:  BIC.mle
# Lookup - 10:  summary.mle
# Lookup - 11:  print.summary.mle
# Lookup - 12:  anova.mle
# Lookup - 13:  print.anova.mle

# Add details/reference to summary method for
#   Wald test

# Lookup - 01
#' Print mle objects
#'
#' @param object a mle object.
#' @export

print.mle = function( object ) {

  if ( !is.null( object$convergence ) ) {
    cat( object$convergence, "\n" )
  }

  if ( !is.null( object$coefficients ) ) {
    K = length( object$coefficients )

    if ( K == 1 )
      cat( paste( K, "estimate:" ), "\n" )
    else
      cat( paste( K, "estimates:" ), "\n" )
    tmp = rbind( object$coefficients )
    rownames( tmp ) = ""
    printCoefmat( tmp )
  }

  cat( "Sum of log-likelihoods:", "\n" )
  cat( round( object$value, 2 ) )
}

# Lookup - 02
#' Extract parameters from mle objects
#'
#' @param object a mle object.
#' @export

coef.mle = function( object ) {
  return( object$coefficients )
}

# Lookup - 03
#' Extract standard errors from mle objects
#'
#' @param object a mle object.
#' @export

sd.mle = function( object ) {
  return( object$SE )
}

# Lookup - 04
#' Extract the variance-covariance matrix from mle objects
#'
#' @param object a mle object.
#' @export

vcov.mle = function( object ) {

  fisher_info = NULL

  if ( !is.null( object$hessian ) ) {

    fisher_info = tryCatch(
      solve(-object$hessian),
      error = function(e) return( NULL ) )
  }

  if ( !is.null( fisher_info ) )
    return( fisher_info )
  else stop( "Unable to estimate variance-covariance matrix" )
}

# Lookup - 05
#' Extract confidence intervals from mle objects
#'
#' @param object a mle object.
#' @export

confint.mle = function( object ) {
  return( object$CI )
}

# Lookup - 06
#' Extract log-likelihood from mle objects
#'
#' @param object a mle object.
#' @export

logLik.mle = function( object ) {

  # Extract sum of the log-likelihoods
  out = object$value

  # Extract number of estimated parameters
  K = length( object$coefficients )

  # Extract number of observations
  n = Compute_N_obs( object )

  # Set degrees of freedom attribute
  attr( out, "df" ) = n - K

  # Set class of output
  class( out ) = "logLik"

  return( out )
}

# Lookup - 07
#' Extract model deviance from mle objects
#'
#' @param object a mle object.
#' @export

deviance.mle = function( object ) {
  return( -2 * object$value )
}

# Lookup - 08
#' Compute AIC for mle objects
#'
#' @param object a mle object.
#' @param finite logical; if true, applies correction for
#'   finite sample size (e.g., Wagenmakers & Farrell, 2004).
#'
#' @details The correction for finite samples is \eqn{2K(K+1)/(n-K-1)}
#'   where K is the number of parameters and n is the number of
#'   observations.
#'
#' @section References:
#' Wagenmakers, E. J., & Farrell, S. (2004). AIC model selection
#'   using Akaike weights. Psychonomic bulletin & review, 11,
#'   192-196.
#'
#' @export

AIC.mle = function( object, finite = F  ) {

  # Extract number of parameters
  K = length( object$coefficients )

  if ( finite ) {
    # Extract number of observations
    n = Compute_N_obs( object )

    # Compute Akaike's information criterion
    # with a finite sample correction
    adjust = ( 2 * K * ( K + 1 ) ) / ( n - K - 1 )
    aic = -2 * object$value + 2 * K + adjust
  } else {
    # Compute Akaike's information criterion
    aic = -2 * object$value + 2 * K
  }

  return( aic )
}

# Lookup - 09
#' Compute BIC for mle objects
#'
#' @param object a mle object.
#' @export

BIC.mle = function( object ) {

  # Extract number of parameters
  K = length( object$coefficients )

  # Extract number of observations
  n = Compute_N_obs( object )

  bic = -2 * object$value + K * log( n )

  return( bic )
}

# Lookup - 10
#' Computes summary measures for mle objects
#'
#' @param object a mle object..
#' @param crit comparison value(s) used in the Wald test.
#' @export

summary.mle = function( object, crit = 0 ) {

  # Details for the sum of the log-likelihoods
  aic = AIC( object )
  bic = BIC( object )
  logLik = object$value

  # Extract number of estimated parameters
  K = length( object$coefficients )

  # Extract information about estimation
  conv = object$convergence
  runTime = object$Time

  # Extract number of observations
  n = Compute_N_obs( object )

  if ( !is.null( object$coefficients ) ) {

    coef = object$coefficients
    est = cbind( Estimates = coef )

    if ( !is.null( object$SE ) ) {

      se =  object$SE
      ci = object$CI

      est = cbind( est,
                   SE = se,
                   ci[1,],
                   ci[2,] )
      colnames( est ) = c(
        "Estimates",
        "SE",
        rownames( ci ) )

      # Carry out Wald tests on parameters

      df_sc = coef - crit
      z = df_sc/se

      p_l = pnorm( z, lower.tail = T )
      p_u = pnorm( z, lower.tail = F )
      pval = pmin( p_l, p_u )

      est = cbind( est, Z = z, p = pval )
    }

  } else est = NULL

  res = list(
    est = est,
    n = n,
    K = K,
    logLik = logLik,
    aic = aic,
    bic = bic,
    conv = conv,
    time = runTime )

  class(res) = "summary.mle"

  return(res)
}

# Lookup - 11
#' Print summary measures for mle objects
#'
#' @param res output from summary of a mle object.
#' @export

print.summary.mle = function( res ) {

  # Convergence and run-time
  if ( res$K > 0 ) {
    string = paste( "Estimation: ", res$conv, ", ",
                    round( res$time[[1]], 2 ), " ",
                    units( res$time ), sep = "" )
  } else {
    string = paste( "Estimation: ",
                    round( res$time[[1]], 2 ), " ",
                    units( res$time ), sep = "" )
  }
  cat( string,  "\n", "\n" )

  if (res$K == 1) nPar = "1 parameter" else
    nPar = paste( res$K, "parameters" )

  cat(
    paste( res$n, "observations," ),
    nPar,
    "\n"
  )

  cat( "\n" )
  cat( "Sum of the log-likelihoods:",
       round( res$logLik, 2 ), "\n" )
  if ( res$K > 0 ) {
    cat(
      paste( "AIC: ", round( res$aic, 2 ), ", ",
             "BIC: ", round( res$bic, 2 ), sep = "" ),
      "\n" )
  } else {
    cat( paste( "Deviance: ", round( res$aic, 2 ),
                sep = "" ), "\n" )
  }

  if ( !is.null( res$est ) ) {

    cat( "\n" )

    if ( ncol( res$est ) > 1 ) {
      printCoefmat( res$est, digits = 2, P.values = T,
                    has.Pvalue = T,
                    signif.stars = F )
    } else {
      printCoefmat( res$est, digits = 2 )
    }

  }

}

# Lookup - 12
#' Table of model comparisons using AIC and BIC
#'
#' @param ... A set of mle objects.
#' @param finite logical; if true, applies the finite sample
#'   correction to AIC calculations.
#' @param modelNames an optional character vector of labels
#'   for the mle objects.
#' @export

anova.mle = function( ..., finite = T, modelNames = NULL ) {

  dots = list( ... )

  # Extract number of models
  M = length( dots )

  # Initialize output
  out = matrix( NA, M, 7 )
  colnames( out ) = c( "N", "K", "Deviance", "AIC",
                       "AIC (w)", "BIC", "BIC (W)" )

  for ( m in 1:M ) {

    out[m,1] = Compute_N_obs( dots[[m]] )
    out[m,2] = length( dots[[m]]$coefficients )
    out[m,3] = deviance( dots[[m]] )
    out[m,4] = AIC( dots[[m]], finite = finite )
    out[m,6] = BIC( dots[[m]] )

  }

  out[,5] = akaike_weights( out[,4] )
  out[,7] = akaike_weights( out[,6] )

  if ( is.null( modelNames ) |
       length( modelNames ) != M ) {
    modelNames = paste( "M", 1:M, sep = "" )
  }
  rownames( out ) = modelNames

  out = list( table = out )
  class( out ) = "anova.mle"

  return( out )
}

# Lookup - 13
#' Print model comparisons for mle objects
#'
#' @param out output from \code{anova} call on a set of mle
#'   objects.
#' @export

print.anova.mle = function( out ) {

  printCoefmat( out[[1]], digits = 2,
                zap.ind = c(5,7) )

}
