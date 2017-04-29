#-----------------------------------#
# Convenience functions for maximum #
# likelihood estimation             #
#-----------------------------------#

# Package development
# library(devtools)
# library(roxygen2)

### TO DO ###
# 1) Build test script
# 2) Optimize param_est function

# Index
# Lookup - 01:  tranPar
# Lookup - 02:  mle_fn_template
# Lookup - 03:  create_mle_fn
# Lookup - 04:  akaike_weights
# Lookup - 05:  param_est'
# Lookup - 06:  MLE*'

# *Needs documentation
# 'Needs examples

# Lookup - 01
#' Transform a vector of parameters
#'
#' A function to transforms or reverse transformations of
#' parameter estimates, which is particularly useful to
#' bound estimates to be positive or restricted to a set
#' range.
#'
#' @param prm a vector of parameters.
#' @param type a vector of the type of transformation to apply
#'   to each parameter, where...
#'   \itemize{
#'     \item{0 = no transformation;}
#'     \item{1 = exponentiation;}
#'     \item{2 = the logistic transform;}
#'     \item{3 = the inverse.}
#'   }
#' @param reverse a logical value; if true, reverses the transformation
#'   specified by \code{type}.
#' @return A vector of transformed or reversed-transformed parameters.
#'
#' @examples
#' prm = c( 1, -1, 0, .5 )
#' prm_star = tran_par( prm, 0:3 )
#' prm_star
#' prm == tran_par( prm_star, 0:3, reverse = T )
#'
#' @export

tran_par = function( prm, type, reverse = F ) {

  if ( !reverse ) {
    prm[ type == 1 ] = exp( prm[ type == 1 ] )
    prm[ type == 2 ] = 1/( 1 + exp( -prm[ type == 2 ] ) )
    prm[ type == 3 ] = 1/prm[ type == 3 ]
  } else {
    prm[ type == 1 ] = log( prm[ type == 1 ] )
    prm[ type == 2 ] = log( prm[ type == 2 ]/(1-prm[ type == 2 ]) )
    prm[ type == 3 ] = 1/prm[ type == 3 ]
  }

  return( prm )
}

# Lookup - 02
#' Template a log-likelihood function
#'
#' Generates a template for a function to compute the log-likelihood
#' or sum of log-likelihoods for data given a vector of parameters.
#' This function can then be passed into the \code{MLE} function.
#'
#' @param version Indicates the type of template to return.
#' @return Prints a template for a function to the console that
#'   may then be copied and modified for the user's needs.
#'
#' @examples
#' mle_fn_template()
#'
#' @export

mle_fn_template = function( version = 0 ) {

  if ( version == 0 ) {
    cat( "mle_fn = function( prm, dat, ", "\n" )
    cat( "                   sum = T, priors = NULL ) {", "\n" )
    cat( "  # Initial setup", "\n" )
    cat( "  ...", "\n" )
    cat( "\n" )
    cat( "  # Calculate the log-likelihoods", "\n" )
    cat( "  ll = f( dat, prm, log = T )", "\n" )
    cat( "  if ( !sum ) return( ll )", "\n" )
    cat( "\n" )
    cat( "  # Sum the log-likelihoods", "\n" )
    cat( "  sll = sum( ll )", "\n" )
    cat( "\n" )
    cat( "  # Incorporate priors (Optional)", "\n" )
    cat( "  if ( !is.null(priors) ) {", "\n" )
    cat( "    sll = sll + ", "\n" )
    cat( "    ...", "\n" )
    cat( "  }", "\n" )
    cat( "\n" )
    cat( "  # Check for NA values", "\n" )
    cat( "  if ( is.na( sll ) ) sll = -Inf", "\n" )
    cat( "\n" )
    cat( "  return( sll )", "\n" )
    cat( "}", "\n" )
  }
}

# Lookup - 03
#' Create a MLE function
#'
#' Given a simple log-likelhood function, produces a new
#' function that can be passed into the \code{MLE} function.
#'
#' @param f a function of the form \code{f(dat,prm)}, where
#'   \code{dat} is the data to be fitted and \code{prm} is
#'   a vector of parameters.
#' @return A new function that takes 4 variables:
#'   \describe{
#'     \item{\code{dat}}{The data to be fitted}
#'     \item{\code{prm}}{A vector of parameters}
#'     \item{\code{sum}}{A logical value for summing the
#'       log-likeihoods}
#'     \item{\code{priors}}{A nuisance variable included for
#'       compatibility}
#'   }
#'
#' @examples
#' f = function(dat,prm) dnorm( dat, prm[1], exp( prm[2] ), log = T )
#' mle_fn = create_mle_fn(f)
#' MLE( rnorm( 100 ), mle_fn, c(0,1) )
#'
#' @export

create_mle_fn = function( f ) {

  mle_fn = function( prm, dat,
                     sum = T, priors = NULL ) {

    ll = f( dat, prm )
    if ( !sum ) return( ll )
    sll = sum( ll )
    if ( is.na( sll ) ) sll = -Inf
    return( sll )

  }

  return( mle_fn )
}

# Lookup - 04
#' Compute Akaike/Schwarz weights
#'
#' Computes Akaike or Schwarz weights (e.g., Wagenmakers & Farrell,
#' 2004) given a vector of information criterion values.
#'
#' @param ic_values a vector of AIC, BIC, or other information
#'   criterion values.
#' @param modelNames an optional vector of labels for each model.
#'
#' @section References:
#' Wagenmakers, E. J., & Farrell, S. (2004). AIC model selection
#'   using Akaike weights. Psychonomic bulletin & review, 11,
#'   192-196.
#'
#' @return A vector of the relative probabilities for each
#'   model.
#'
#' @examples
#' data( ToothGrowth ) # Load in data
#' # Fit 2 regression models
#' m1 = lm( len ~ dose, data = ToothGrowth )
#' m2 = lm( len ~ dose + supp, data = ToothGrowth )
#' # Create vector of AIC values
#' ic = c( AIC(m1), AIC(m2) )
#' akaike_weights( ic )
#' # Include meaningful labels
#' modelNames = c( 'Dose', 'Dose+Supp' )
#' akaike_weights( ic, modelNames )
#'
#' @export

akaike_weights = function( ic_values,
                           modelNames = NULL ) {

  delta_ic = ic_values - min( ic_values )
  rel_like = exp( -.5 * delta_ic )
  w_ic = rel_like / sum( rel_like )

  if ( is.null( modelNames ) |
       length( modelNames ) != length( ic_values ) ) {
    modelNames = paste( "M", 1:length( ic_values ), sep = "" )
  }
  names( w_ic ) = modelNames

  return( w_ic )
}

# Lookup - 05
#' Compute parameter matrix
#'
#' A function to calculate the weighted sum of a design matrix
#' and a set of coefficients, producing a matrix of parameters
#' by observations. Fixed coefficient values can also be specified.
#'
#' @param X a design matrix.
#' @param coef a vector of C coefficients.
#' @param fixed a vector of F values for the fixed coefficients.
#' @param index a matrix with two columns, the first giving the
#'  row positions for the free coefficients (followed by the
#'  row positions for any fixed values), the second giving the
#'  column positions for the free coefficients (followed by
#'  the positions for the fixed values). The first column of
#'  the final row of the matrix indicates the total number of
#'  desired parameters so that the parameter matrix can be
#'  created.
#' @param parSel a vector giving the indices mapping the values
#'   in the coefficient vector to the rows of the 'index' matrix,
#'   thereby allowing different conditions to be constrained
#'   to have the same free coefficient.
#' @details Given N observations and P desired parameters, the goal
#'   is to produce a P x N matrix given a set of V covariates, C
#'   coefficients, and F fixed values. To do so, a P x V parameter
#'   matrix M is specified, and the 'index' matrix along with the
#'   \code{parSel} vector are used to fill the positions of the
#'   P x V matrix. Fixed values from the \code{fixed} vector are
#'   additionally included in the matrix M. To produce the
#'   desired P x N output matrix, the P x V matrix M is multiplied by
#'   the V x N design matrix X.
#' @return A P x N matrix giving the set of parameter values for
#'   each of the N observations.
#'
#' @export

param_est = function( X, coef, fixed, index,
                      parSel ) {

  N = ncol(X)
  V = nrow(X)
  C = length( parSel )
  L = dim(index)[1]
  P = index[ L, 1 ]
  finish = L - 1 - C

  pm = matrix( 0.0, P, V )

  pm[ matrix( index[ 1:C, ], C, 2 ) ] = coef[parSel]

  if (finish > 0) {
    pm[ matrix( index[ 1:finish + C, ], finish, 2 ) ] = fixed
  }

  return( pm %*% X )
}

# Lookup - 06
#' A convenience function to carry out maximum likelihood estimation
#' using either \code{optim} or \code{nlm}.
#'
#' @param dat the data to be fitted; can be in a variety of forms
#'   as long as it is compatible with the \code{mle_fn} function
#'   that is passed in.
#' @param mle_fn a function to compute the sum of the log-likelihoods
#'   for a vector of parameters \code{prm}, the data \code{dat}, and
#'   control parameters \code{sum} and \code{priors}.
#' @param start either a vector of initial starting values for the
#'   parameters, or a function (which takes no parameters) used to
#'   generate a dispersed set of starting values.
#' @param gr_fn an optional function to compute the gradients for the
#'   parameters being estimated. See \code{\link{optim}}.
#' @param method the optimization algorithm to use. See
#'   \code{\link{optim}}. If \code{NULL}, a step-by-step approach
#'   is used, where viable estimates and standard errors are sought
#'   by first trying the 'BFGS' method, then the 'BFGS' method using
#'   initial estimates from a 'Nelder-Mead' approach, then the 'nlm'
#'   approach, and finally the 'Nelder-Mead' algorithm.
#' @param control a named list where options controlling the
#'   optimization can be specified. See \code{\link{optim}} and
#'   \code{\link{nlm}}. The lower and upper limits for the 'Brent'
#'   method can also be specified in this list.
#' @param hessian a logical value; if true, the hessian matrix is
#'   estimated, allowing the derivation of standard errors and
#'   confidence intervals.
#' @param alpha the width of the interval for the confidence intervals
#'   around parameters. Defaults to +/- two standard deviations.
#' @param emStop the number of iterations used when attempting to
#'   generate viable starting values.
#' @param nRep the number of times to repeat estimation using
#'   dispersed starting values. Only applicable when a function to
#'   generate starting values is provided.
#' @param parNames an optional character vector giving a set of labels
#'   for the parameters being estimated.
#' @return Forthcoming.
#'
#' @export

MLE = function(dat, mle_fn, start,
               priors = NULL,
               gr_fn = NULL,
               method = NULL,
               control = list(),
               hessian = T,
               alpha = .9544997,
               emStop = 20,
               nRep = 1,
               parNames = NULL ) {

  # Track run-time
  start_time = Sys.time()

  # Initialize variables
  track_value = rep( -Inf, nRep )
  prev_value = -Inf
  out = list()

  # If starting values are not generated by
  # a function, set iterations to 1
  if ( is.vector( start ) | is.null( start ) ) nRep = 1

  # Loop over iterations
  for ( nr in 1:nRep ) {
    output = tryCatch(
      mleWrapper( dat, mle_fn, start,
                  priors = priors,
                  gr_fn = gr_fn,
                  method = method,
                  control = control,
                  hessian = hessian,
                  alpha = alpha,
                  emStop = emStop ),
      error = function(e) {
        message(e)
        return( NULL )
      } )

    if ( !is.null( output ) ) {

      track_value[nr] = output$value

      if ( output$value > prev_value ) {
        out = output
        prev_value = output$value
      }
    }
  }

  # Check in case all estimation attempts failed
  if ( sum( track_value == -Inf ) == nRep )
    stop( 'All estimation attempts failed', call. = F )

  # Track previous sums of log-likelihoods
  # (Tracks possibility of local maxima)
  out$track_value = track_value

  # Include maximum likelihood function in output
  out$mle_fn = mle_fn

  # Include prior vlaues
  out$priors = priors

  # Determine number of coefficients
  K = length( out$coefficients )

  if ( K > 0 ) {
    # If no parameter names are specified or
    # if names given do not match number of
    # parameters

    if ( is.null( parNames ) |
         length( parNames ) != K ) {
      parNames = paste( 'V', 1:K, sep = '' )
    }

    # Add variable names to output
    names( out$coefficients ) = parNames
    names( out$start ) = parNames

    if ( !is.null( out$hessian ) ) {
      colnames( out$hessian ) = parNames
      rownames( out$hessian ) = parNames
    }
    if ( !is.null( out$SE ) ) {
      names( out$SE ) = parNames
    }
    if ( !is.null( out$CI ) ) {
      colnames( out$CI ) = parNames
      rownames( out$CI ) = c(
        paste( round( 100 * (1-alpha)/2 ), '%', sep = '' ),
        paste( round( 100 * ( (1-alpha)/2 + alpha ) ), '%', sep = '' )
      )
    }
  }

  run_time = Sys.time() - start_time
  out$Time = run_time

  # Create new class for output
  class(out) = append( class(out), "mle" )

  return( out )
}
