#--------------------#
# Internal functions #
#--------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01:  checkLik
# Lookup - 02:  SE_and_CI
# Lookup - 03:  one_param
# Lookup - 04:  multi_param
# Lookup - 05:  mle_steps
# Lookup - 06:  convergenceChecker
# Lookup - 07:  mleWrapper*
# Lookup - 08:  Compute_N_obs

# *Needs documentation

# Lookup - 01
checkLik = function( start, dat, mle_fn, priors, emStop ) {
  # Purpose:
  # A function to generate starting values that produce
  # valid sums of log-likelihoods.
  # Arguments:
  # start  - A vector of starting values or a function that
  #          generates dispersed starting values
  # dat    - The data to be fitted
  # mle_fn - A function to compute the (possibly penalized)
  #          sum of the log-likelihoods
  # priors - An optional input for priors on the parameters
  # emStop - The number of iterations to attempt to find valid
  #          outputs before stopping
  # Returns:
  # A vector of starting values for parameter estimation.

  # 1) Check for errors
  # 2) Check if feasible starting

  # Fixed starting values
  if ( is.null( start ) | is.vector( start ) ) {
    start_val = start
    emStop = 1
  }

  # Initialize variables
  value = -Inf

  for ( em in 1:emStop ) {

    if ( !is.null( start ) & !is.vector( start ) ) {
      # Generate dispersed starting values
      start_val = start()
    }

    # Compute sum of the log-likeihoods
    chk = tryCatch(
      # Attempt to run function
      mle_fn( start_val, dat, sum = T, priors = priors ),
      # If there is an error
      error = function(e) {
        return( e )
      }
    )

    if ( is.numeric( chk ) ) value = chk

    # Stop the loop if starting values work
    if ( value != -Inf ) break;

  }

  if ( value == -Inf & is.numeric( chk ) & em == emStop )
    stop("No feasible starting values were found \n", call. = F )
  if ( value == -Inf & !is.numeric( chk ) & em == emStop )
    stop( paste( "Error with likelihood;", chk[1], "\n" ), call. = F )

  return( start_val )
}

# Lookup - 02
SE_and_CI = function( par, hessian, alpha ) {
  # Purpose:
  # Estimates the standard error and confidence intervals
  # for a set of parameters based on the hessian matrix
  # returned by the 'optim' or 'nlm' functions.
  # Arguments:
  # par     - A vector of parameter estimates
  # hessian - The estimated hessian matrix
  # alpha   - The width of the confidence interval
  # Returns:
  # A list with the standard errors and a matrix whose
  # rows are the lower and upper limits for the confidence
  # interval; if there were errors or inadmissable values,
  # returns NULL.

  # Initialize output
  out = NULL

  # Attempt to calculate the Fisher information matrix
  fisher_info = tryCatch(
    solve(-hessian),
    error = function(e) return( NULL ) )

  # If matrix is invertable
  if ( !is.null( fisher_info ) ) {

    # Calculate standard error
    prop_sigma = sqrt(diag(fisher_info))

    # Calculate critical value
    if ( !is.numeric( alpha ) |
         ( alpha > 1 | alpha < 0 ) ) alpha = .9544997
    crt = abs( qnorm( (1-alpha)/2 ) )

    # Calculate confidence intervals
    CI = matrix( NA, 2, length( par ) )
    CI[1,] = par - crt*prop_sigma
    CI[2,] = par + crt*prop_sigma

    # If no NA values
    if ( sum( is.na( prop_sigma ) ) == 0 &
         sum( is.na( CI ) ) == 0 ) {
      out = list( SE = prop_sigma,
                  CI = CI )
    }

  }

  return( out )
}

# Lookup - 03
one_param = function( dat, mle_fn, startVal, priors,
                      gr_fn, control, hessian, alpha ) {
  # Purpose:
  # A function to carry out maximum likelihood optimization
  # for a single parameter.
  # Arguments:
  # dat      - The data to be fitted
  # mle_fn   - A function to compute the (possibly penalized)
  #            sum of the log-likelihoods
  # startVal - A vector of starting values to be passed into mle_fn
  # priors   - Optional input for priors on the parameters
  # gr_fn    - An optional function to compute the gradients for the
  #            parameters
  # control  - A named list of values for additional control of the
  #            'optim' function
  # hessian  - A logical value; if true, the hessian matrix is
  #            estimated and returned
  # alpha    - The width for the confidence intervals
  # Returns:
  # The output from the 'optim' function.

  control$fnscale = -1 # Set optimization to maximize

  if ( is.null( control$lower ) ) {
    lower = -100
  } else {
    lower = control$lower; control$lower = NULL
  }

  if ( is.null( control$upper ) ) {
    upper = 100
  } else {
    upper = control$upper; control$upper = NULL
  }

  results = suppressWarnings(
    optim( startVal, mle_fn,
           dat = dat,
           priors = priors,
           sum = T,
           gr = gr_fn,
           method = "Brent",
           lower = lower,
           upper = upper,
           control = control,
           hessian = hessian )
  )

  # Store information about optimization routine
  results$algorithm = list(
    type = 'Brent (optim)',
    counts = results$counts,
    message = results$message )

  # Compute standard errors and confidence intervals
  results$SE = NULL; results$CI = NULL
  if ( hessian ) {
    tmp = SE_and_CI( results$par,
                     results$hessian,
                     alpha )
    if (!is.null(tmp)) {
      results$SE = tmp$SE
      results$CI = tmp$CI
    }
  }

  return( results )
}

# Lookup - 04
multi_param = function( dat, mle_fn, startVal, priors,
                        gr_fn, control, hessian, alpha,
                        method ) {
  # Purpose:
  # A function to carry out maximum likelihood optimization
  # for multiple parameters.
  # Arguments:
  # dat      - The data to be fitted
  # mle_fn   - A function to compute the (possibly penalized)
  #            sum of the log-likelihoods
  # startVal - A vector of starting values to be passed into mle_fn
  # priors   - Optional input for priors on the parameters
  # gr_fn    - An optional function to compute the gradients for the
  #            parameters
  # control  - A named list of values for additional control of the
  #            'optim' function
  # hessian  - A logical value; if true, the hessian matrix is
  #            estimated and returned
  # alpha    - The width for the confidence intervals
  # method   - The optimization algorithm to use
  # Returns:
  # The output from the 'optim' function or output from the 'nlm'
  # function revised to match the 'optim' output.

  if ( method != 'nlm' ) {
    control$fnscale = -1 # Set optimization to maximize

    if ( !is.null( control$lower ) ) {
      lower = control$lower; control$lower = NULL
    } else lower = -Inf

    if ( !is.null( control$upper ) ) {
      upper = control$upper; control$upper = NULL
    } else upper = Inf

    results = suppressWarnings(
      optim( startVal, mle_fn,
             dat = dat,
             priors = priors,
             sum = T,
             gr = gr_fn,
             lower = lower,
             upper = upper,
             method = method,
             control = control,
             hessian = hessian )
    )

    # Store information about algorithm
    results$algorithm = list(
      type = paste( method, '(optim)' ),
      counts = results$counts,
      message = results$message
    )

    # Compute standard errors and confidence intervals
    results$SE = NULL; results$CI = NULL
    if ( hessian ) {
      tmp = SE_and_CI( results$par,
                       results$hessian,
                       alpha )
      if (!is.null(tmp)) {
        results$SE = tmp$SE
        results$CI = tmp$CI
      }
    }

  } else {

    # Change to minimization problem
    mle_fn_min = function( prm, dat, priors )
      return( -mle_fn( prm, dat, sum = T, priors = priors ) )

    results = suppressWarnings( nlm( mle_fn_min, startVal,
                                     dat = dat, priors = priors,
                                     hessian = hessian ) )

    # Convert results to match output of 'optim'
    results$convergence = 0
    if ( results$code > 2 )
      results$convergence = 1
    results$code = NULL
    results$par = results$estimate
    results$estimate = NULL
    results$value = -results$minimum
    results$minimum = NULL
    results$counts = list(
      counts = results$iterations,
      gradients = results$gradient
    )
    results$gradient = NULL
    if ( hessian ) {
      K = dim( results$hessian )[1]
      for ( k in 1:K )
        results$hessian[k,k] =
          -results$hessian[k,k]
    }

    # Store information about algorithm
    results$algorithm = list(
      type = 'nlm',
      counts = results$counts,
      message = NULL
    )

    # Compute standard errors and confidence intervals
    results$SE = NULL; results$CI = NULL
    if ( hessian ) {
      tmp = SE_and_CI( results$par,
                       results$hessian,
                       alpha )
      if (!is.null(tmp)) {
        results$SE = tmp$SE
        results$CI = tmp$CI
      }
    }

  }

  return( results )
}

# Lookup - 05
mle_steps = function( dat, mle_fn, startVal, priors,
                      gr_fn, control, hessian, alpha ) {
  # Purpose:
  # A function to carry out maximum likelihood optimization
  # for multiple parameters, using a step-by-step process in which
  # several different optimization algorithms are attempted until
  # viable estiamtes are obtained.
  # Arguments:
  # dat      - The data to be fitted
  # mle_fn   - A function to compute the (possibly penalized)
  #            sum of the log-likelihoods
  # startVal - A vector of starting values to be passed into mle_fn
  # priors   - Optional input for priors on the parameters
  # gr_fn    - An optional function to compute the gradients for the
  #            parameters
  # control  - A named list of values for additional control of the
  #            'optim' function
  # hessian  - A logical value; if true, the hessian matrix is
  #            estimated and returned
  # alpha    - The width for the confidence intervals
  # Returns:
  # The output from the 'optim' function or output from the 'nlm'
  # function revised to match the 'optim' output.

  # Check for convergence and viable estimates
  convCheck = function( results ) {

    # Initialize output
    out = FALSE

    # If output was returned
    if ( !is.null( results ) ) {

      # If the model converged
      if ( results$convergence == 0 ) {

        # If standard errors and confidence
        # intervals were supposed to be estimated
        if ( hessian ) {

          if ( !is.null( results$SE ) &
               !is.null( results$CI ) ) out = TRUE

        } else {
          out = TRUE
        }

      }
    }

    return( out )
  }

  # Step 1) Use method BFGS (optim)
  results = tryCatch(
    multi_param( dat, mle_fn, startVal,
                 priors, gr_fn, control, hessian, alpha,
                 method = 'BFGS' ),
    error = function(e) return(NULL)
  )

  # Check for convergence
  if ( convCheck( results ) ) return( results )

  # Otherwise...
  # Step 2) Try Nelder-Mead, then try BFGS (optim)
  results = tryCatch(
    {
      # First, obtain estimates via Nelder-mead
      if ( !is.null( control$maxit ) ) {
        orig_maxit = control$maxit
        control$maxit = 500
      }
      results = multi_param( dat, mle_fn, startVal,
                             priors, gr_fn, control, hessian,
                             method = 'Nelder-mead' )

      # Then, try BFGS using estimates from Nelder-Mead
      if ( exists( "orig_maxit" ) ) {
        control$maxit = orig_maxit
      }
      results = multi_param( dat, mle_fn, results$par,
                             priors, gr_fn, control, hessian,
                             method = 'BFGS' )
    },
    error = function(e) return(NULL)
  )

  # Check for convergence
  if ( convCheck( results ) ) return( results )

  # Otherwise...
  # Step 3) Try 'nlm'
  results = tryCatch(
    results = multi_param( dat, mle_fn, results$par,
                           priors, gr_fn, control, hessian,
                           method = 'nlm' ),
    error = function(e) return(NULL)
  )

  # Check for convergence
  if ( convCheck( results ) ) return( results )

  # Otherwise
  # Step 4) Try Nelder-Mead (optim)
  results = tryCatch(
    results = multi_param( dat, mle_fn, startVal,
                           priors, gr_fn, control, hessian,
                           method = 'Nelder-mead' ),
    error = function(e) return(NULL)
  )

  # Check for convergence
  if ( convCheck( results ) )
    return( results )
  else
    stop( 'Estimation failed' )
}

# Lookup - 06
convergenceChecker = function( code ) {
  # Purpose:
  # Converts the convergence codes from 'optim' into
  # interpretable messages.
  # Arguments:
  # code - the convergence code from 'optim' output.
  # Returns:
  # A character string.

  string = "Unknown error"

  if ( code == -1 ) string = "Fixed parameters"
  if ( code == 0 ) string = "Successful convergence"
  if ( code == 1 ) string = "Maximum number of iterations reached"
  if ( code == 10 ) string = "Degeneracy of Nelder-Mead simplex"
  if ( code == 51 ) string = "Warning for L-BFGS-B method"
  if ( code == 52 ) string = "Error for L-BFGS-B method"

  return( string )
}

# Lookup - 07
mleWrapper = function( dat, mle_fn, start,
                       priors = NULL,
                       gr_fn = NULL,
                       method = NULL,
                       control = list(),
                       hessian = T,
                       alpha = .9544997,
                       emStop = 20 ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  # Initialize output
  output = list(
    coefficients = NULL,
    convergence = NULL,
    start = NULL,
    algorithm = NULL,
    value = NULL,
    hessian = NULL,
    SE = NULL,
    CI = NULL,
    data = dat
  )

  # Define error message
  error_message = function(e) {
    return( NULL )
  }

  # Check starting values
  startVal = checkLik( start, dat, mle_fn, priors, emStop )
  output$start = startVal

  # If parameters are fixed
  if ( is.null( startVal ) ) {

    results = output
    results$value = tryCatch(
      mle_fn( start, dat, sum = T, priors = priors ),
      error = error_message )
    results$convergence = -1
    results$algorithm = "Fixed parameters"

  }

  # If estimating only a single parameter
  if ( length( startVal ) == 1 ) {
    results = tryCatch(
      one_param( dat, mle_fn, startVal, priors, gr_fn,
                 control, hessian, alpha ),
      error = error_message )
  }

  # If estimating multiple parameters
  if ( length( startVal ) > 1 ) {

    # If no method is provided
    if ( is.null( method ) ) {
      results = tryCatch(
        mle_steps( dat, mle_fn, startVal, priors,
                   gr_fn, control, hessian, alpha ),
        error = error_message
      )

    } else {
      results = tryCatch(
        multi_param( dat, mle_fn, startVal, priors,
                     gr_fn, control, hessian, alpha, method ),
        error = error_message
      )
    }

  }

  # Create final output
  if ( !is.null( results ) ) {

    # Extract parameter estimates
    output$coefficients = results$par

    # Store information about convergence
    output$convergence = convergenceChecker(
      results$convergence )

    # Store information about algorithm
    output$algorithm = results$algorithm

    # Store sum of the log-likelihoods
    output$value = results$value

    # Extract standard errors
    output$SE = results$SE

    # Extract confidence intervals
    output$CI = results$CI

  } else {
    stop()
  }

  return( output )
}

# Lookup - 08
Compute_N_obs = function( object ) {
  # Purpose:
  # A function to extract the number of observations that were
  # fitted from a mle object.
  # Arguments:
  # object - a mle object (output from the 'MLE' function)
  # Returns:
  # The number of observations that were fitted.

  n = length(
    object$mle_fn(
      object$coefficients,
      object$dat,
      sum = F, priors = object$priors ) )

  return( n )
}
