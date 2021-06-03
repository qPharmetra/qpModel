globalVariables(c(
'A2periph', 'Aabs', 'Acentral', 'Aperiph', 'CL', 'DV', 'IPRED', 'LF1', 'MTT', 'NTR', 'Q', 'Q2', 'V', 'V2', 'V3',
'abs0', 'add', 'amax', 'cl', 'elim', 'etaCL', 'etaQ', 'etaQ2', 'etaV', 'etaV2', 'etaV3', 'etaamax', 'etacl',
'etaf1', 'etak0', 'etak12', 'etak13', 'etak21', 'etak31', 'etaka', 'etaka50', 'etake', 'etakm50',
'etakmax', 'etamtt', 'etantr', 'inpt', 'k0', 'k12', 'k13', 'k21', 'k31', 'ka', 'ka50', 'ke', 'km50', 'kmax',
'prop', 'str_replace', 'tvCL', 'tvQ', 'tvQ2', 'tvV', 'tvV2', 'tvV3', 'tvamax', 'tvcl', 'tvf1', 'tvk0',
'tvk12', 'tvk13', 'tvk21', 'tvk31', 'tvka', 'tvka50', 'tvke', 'tvkm50', 'tvkmax', 'tvmtt', 'tvntr'
))

#' Create a qpModel

#' Creates an object of class 'qpmodel', convertible
#' to NONMEM code or mrgsolve code.
#'
#' @param ode ordinary differential equations
#' @param algebraic algebraic equations
#' @param param parameters
#' @param omega second level random effects
#' @param sigma first level random effects
#' @param theta fixed effects
#' @param observe items to include in output
#' @param output other output
#' @param global global equations
#' @param custom custom equations
#' @export
#' @family qpmodel
#' @examples
#' qpModel()
#'
qpModel <- function(
  ode = eqns(),
  algebraic = eqns(),
  param = eqns(),
  omega = eqns(),
  sigma = eqns(),
  theta = eqns(),
  observe = eqns(),
  output = eqns(),
  global = eqns(),
  custom = eqns()
){
  x <- list(
    ode = ode,
    algebraic = algebraic,
    param = param,
    omega = omega,
    sigma = sigma,
    theta = theta,
    observe = observe,
    output = output,
    global = global,
    custom = custom
  )

  class(x) <- 'qpmodel'
  return(x)
}

#' Generate Equations
#'
#' Generates equations. See vignettes.
#'
#' @param ... unquoted arguments
#' @export
#' @importFrom rlang enexprs
#' @family eqns
#' @examples
#' eqns()
eqns <- function(...){
  exprns <- enexprs(...)
  structure(list(calls = exprns), class = 'eqns')
}
#' Make an Equation
#'
#' Makes a single equation of class 'eqn'.
#'
#' @param arg instance of class 'call'
#' @family eqn

eqn.make <- function(arg) {
  # arg is a class::call

  # this creates a string from the call for operations
  instring = deparse(arg)
  structure(list(
    call = arg,
    instring = instring,
    vars = all.vars(arg)
  ),
  class = 'eqn')
}

#' Convert Equation to a Call
#'
#' Converts an equation to a call.
#'
#' @param arg an equation
#' @family eqn
#' @export
#'
eqn_ <- function(arg) {
  # arg is a expression

  # this converts arg into a call
  arg$new <- arg
  arg <- arg[[1]]
  eqn.make(arg)
  #structure(list(call = arg),class='eqn')
}

#' Convert Equations to Calls
#'
#' Converts equations to calls.
#' @param ... equations
#' @family eqns
eqns_ <- function(...) {
  # arg is a expression
  exprns = enexprs(...)

  # this converts arg into a call
  exprns$new <- exprns
  exprns <- exprns[length(match.call) - 1]
  #eqn.make(exprns)
  structure(list(calls = exprns), class = 'eqns')
}

#' Add to Equation
#'
#' Adds to equations.
#' @export
#' @family operators
#' @param eq1 equation 1
#' @param eq2 equation 2
`+.eqn` <- function(eq1, eq2) {
  neweq = parse(text = paste('(', eq1$instring, ')+(', eq2$instring, ')'))
  eqn_(neweq)

}

#' Subtract from Equation
#'
#' Subtracts from equations
#' @export
#' @family operators
#' @param eq1 equation 1
#' @param eq2 equation 2
`-.eqn` <- function(eq1, eq2) {
  neweq = parse(text = paste('(', eq1$instring, ')-(', eq2$instring, ')'))
  eqn_(neweq)

}

#' Multiply Equation
#'
#' Multiplies equation.
#' @export
#' @family operators
#' @param eq1 equation 1
#' @param eq2 equation 2
`*.eqn` <- function(eq1, eq2) {
  neweq = parse(text = paste('(', eq1$instring, ')*(', eq2$instring, ')'))
  eqn_(neweq)

}

#' Divide Equation
#'
#' Divides equation.
#' @export
#' @family operators
#' @param eq1 equation 1
#' @param eq2 equation 2
`/.eqn` <- function(eq1, eq2) {
  neweq <- parse(text = paste('(', eq1$instring, ')/(', eq2$instring, ')'))
  eqn_(neweq)


}

#' Exponentiate Equation
#'
#' Exponentiates equation.
#' @export
#' @family operators
#' @param eq1 equation 1
#' @param eq2 equation 2

`^.eqn` <- function(eq1, eq2) {
  neweq <- parse(text = paste('(', eq1$instring, ')^(', eq2$instring, ')'))
  eqn_(neweq)

}

#' Add to Equations
#'
#' Adds to equations.
#' @export
#' @family operators
#' @param eqns1 equation 1
#' @param eqns2 equation 2
`+.eqns` <- function(eqns1, eqns2) {
  # set up new eqns object
  neweqns <- eqns()

  # find calls in both eqn structures
  interx <- intersect(names(eqns1$calls), names(eqns2$calls))

  # find location of duplicates in eqns1
  locate <- names(eqns1$calls) %in% interx

  # create addition of two objects
  neweqns$calls <- c(eqns1$calls[!locate], eqns2$calls)

  structure(list(calls = neweqns$calls), class = 'eqns')

}



#' Multiply Equations
#'
#' Multiplies equations.
#' @export
#' @family operators
#' @param eqns1 equation 1
#' @param eqns2 equation 2
`*.eqns` <- function(eqns1, eqns2) {
  # set up new eqns object
  neweqns <- eqns()

  # create addition of two objects
  neweqns$calls <- c(eqns1$calls, eqns2$calls)

  structure(list(calls = neweqns$calls), class = 'eqns')

}

#' Add 'qpmodel'
#'
#' Adds 'qpmodel'.
#' @export
#' @family operators
#' @param qp1 qpmodel 1
#' @param qp2 qpmodel 2
`+.qpmodel` <- function(qp1, qp2) {
  #odesnew <- append(qp1$ode$calls,qp2$ode$calls)
  #algebraicnew <- append(qp1$algebraic$calls,qp2$algebraic$calls)
  #paramsnew <- append(qp1$param$calls,qp2$param$calls)
  #thetanew <- append(qp1$theta$calls,qp2$theta$calls)
  #omeganew <- append(qp1$omega$calls,qp2$omega$calls)
  #sigmanew <- append(qp1$sigma$calls,qp2$sigma$calls)
  qpnew <- qpModel()

  qpnew$ode <- qp1$ode + qp2$ode
  qpnew$algebraic <- qp1$algebraic + qp2$algebraic
  qpnew$param <- qp1$param + qp2$param
  qpnew$theta <- qp1$theta + qp2$theta
  qpnew$omega <- qp1$omega + qp2$omega
  qpnew$sigma <- qp1$sigma + qp2$sigma
  qpnew$observe <- qp1$observe + qp2$observe
  qpnew$output <- qp1$output + qp2$output
  qpnew$global <- qp1$global + qp2$global
  qpnew$custom <- qp1$custom + qp2$custom

  #odesnew <- eqn_(odesnew)

  #qp1$ode <- qp1$ode + qp2$ode
  #qp1$algebraic <- qp1$algebraic + qp2$algebraic
  #qp1$param <- qp1$param + qp2$param
  #qp1$theta <- qp1$theta + qp2$theta
  #qp1$omega <- qp1$omega + qp2$omega
  #qp1$sigma <- qp1$sigma + qp2$omega

  qpModel(
    qpnew$ode,
    qpnew$algebraic,
    qpnew$param,
    qpnew$omega,
    qpnew$sigma,
    qpnew$theta,
    qpnew$observe,
    qpnew$output,
    qpnew$global,
    qpnew$custom
  )
}

#' Muliply 'qpmodel'
#'
#' Multiplies 'qpmodel'.
#' @export
#' @family operators
#' @param qp1 qpmodel 1
#' @param qp2 qpmodel 2

`*.qpmodel` <- function(qp1, qp2) {
  #odesnew <- append(qp1$ode$calls,qp2$ode$calls)
  #algebraicnew <- append(qp1$algebraic$calls,qp2$algebraic$calls)
  #paramsnew <- append(qp1$param$calls,qp2$param$calls)
  #thetanew <- append(qp1$theta$calls,qp2$theta$calls)
  #omeganew <- append(qp1$omega$calls,qp2$omega$calls)
  #sigmanew <- append(qp1$sigma$calls,qp2$sigma$calls)
  qpnew <- qpModel()

  qpnew$ode <- qp1$ode * qp2$ode
  qpnew$algebraic <- qp1$algebraic * qp2$algebraic
  qpnew$param <- qp1$param * qp2$param
  qpnew$theta <- qp1$theta * qp2$theta
  qpnew$omega <- qp1$omega * qp2$omega
  qpnew$sigma <- qp1$sigma * qp2$sigma
  qpnew$observe <- qp1$observe * qp2$observe
  qpnew$output <- qp1$output * qp2$output
  qpnew$global <- qp1$global * qp2$global
  qpnew$custom <- qp1$custom * qp2$custom

  #odesnew <- eqn_(odesnew)

  #qp1$ode <- qp1$ode + qp2$ode
  #qp1$algebraic <- qp1$algebraic + qp2$algebraic
  #qp1$param <- qp1$param + qp2$param
  #qp1$theta <- qp1$theta + qp2$theta
  #qp1$omega <- qp1$omega + qp2$omega
  #qp1$sigma <- qp1$sigma + qp2$omega

  qpModel(
    qpnew$ode,
    qpnew$algebraic,
    qpnew$param,
    qpnew$omega,
    qpnew$sigma,
    qpnew$theta,
    qpnew$observe,
    qpnew$output,
    qpnew$global,
    qpnew$custom
  )
}

#' View 'qpmodel'
#'
#' Views 'qpmodel'.
#' @export
#' @param qpmodel an object of class 'qpmodel'
#' @family qpmodel
modelView <- function(qpmodel = qpModel()) {
  stopifnot(inherits(qpmodel, 'qpmodel'))
  # write out ode
  #print('ODEs:')
  #print(paste(names(qpmodel$ode$calls),' = ',qpmodel$ode$calls))
  ode <- data.frame('ODEs' = paste(names(qpmodel$ode$calls), ' = ', qpmodel$ode$calls))

  # write out relationships
  #print('Algebraic eqns:')
  #print(paste(names(qpmodel$algebraic$calls),' = ',qpmodel$algebraic$calls ))
  algebra <- data.frame('Algebraic' = paste(
    names(qpmodel$algebraic$calls),
    ' = ',
    qpmodel$algebraic$calls
  ))

  # write out parameters
  #print('Parameters:')
  #print(paste(names(qpmodel$param$calls),' = ',qpmodel$param$calls ))
  param <- data.frame('Parameters' = paste(names(qpmodel$param$calls), ' = ', qpmodel$param$calls))

  # write out fixed effect values
  #print('Fixed Effect values:')
  #print(paste(names(qpmodel$theta$calls),' = ',qpmodel$theta$calls ))
  theta <- data.frame('Theta' = paste(names(qpmodel$theta$calls), ' = ', qpmodel$theta$calls))

  # write out omega values
  #print('Omega:')
  #print(paste(names(qpmodel$omega$calls),' = ',qpmodel$omega$calls ))
  omega <- data.frame('Omega' = paste(names(qpmodel$omega$calls), ' = ', qpmodel$omega$calls))

  # write out observation info
  #print('Observations:')
  #print(paste(names(qpmodel$observe$calls),' = ',qpmodel$observe$calls ))
  observe <- data.frame('Observe' = paste(names(qpmodel$observe$calls), ' = ', qpmodel$observe$calls))

  # write out error values
  #print('Residual Error:')
  #print(paste(names(qpmodel$sigma$calls),' = ',qpmodel$sigma$calls ))
  sigma <- data.frame('Sigma' = paste(names(qpmodel$sigma$calls), ' = ', qpmodel$sigma$calls))

  # write out output variables
  output <- data.frame('Output' = paste(qpmodel$output$calls))

  preview <- list(
    'ODEs' = ode,
    'Algebraic' = algebra,
    'Params' = param,
    'Theta' = theta,
    'Omega' = omega,
    'Observe' = observe,
    'Sigma' = sigma,
    'Output' = output
  )

  return(preview)
}

#' Coerce to mrgsolve
#'
#' Coerces model to mrgsolve format
#'
#' @param qpmodel class 'qpmodel'
#' @param filename output filename
#' @export
#'
modelCreate_R <- function(
  qpmodel = qpModel(),
  filename = 'model.txt'
){
  ## Create model text file
  filemodel <- file(filename)

  ## start writing model code for mrgsolve
  #writeLines(c('modelfile =  '), filemodel)

  ## create parameter list
  finalparams <- paste(names(qpmodel$theta$calls), '=', qpmodel$theta$calls)

  ## create compartment list
  finalcmts <- names(qpmodel$ode$calls)

  ## create main
  finalmain <- paste(
    'double',
    names(qpmodel$param$calls),
    '=',
    qpmodel$param$calls,
    ';'
  )

  ## create ODEs
  ode <- paste(
    'dxdt_',
     names(qpmodel$ode$calls),
     ' = ',
     qpmodel$ode$calls,
     ';',
     sep = ''
  )

  ## create Algebraic eqns
  algebra <- character(0)

  if(length(qpmodel$algebraic$calls)){
    algebra <- paste(
      'double',
      names(qpmodel$algebraic$calls),
      '=',
      qpmodel$algebraic$calls,
      ';'
    )
  }

  ## create SIGMA
  sigma <- paste(qpmodel$sigma$calls)

  ## create OMEGA
  omega <- paste(qpmodel$omega$calls)

  ## create observation equations
  if (is.null(qpmodel$observe$call)){
    table <- NULL
  } else{
    table <- paste(
      'double',
      names(qpmodel$observe$calls),
      ' = ',
      qpmodel$observe$calls,
      ';'
    )
  }

  ## create table output
  capture <- paste(qpmodel$output$calls)

  ### TCAM components

  ## MAIN & ODE
  if (is.null(qpmodel$custom$calls)) {
    tcammain <- NULL
    tcamode <- NULL
  } else{
    tcammain <- qpmodel$custom$calls$r.main
    tcamode <- qpmodel$custom$calls$r.ode
  }

  ## GLOBAL
  if (is.null(qpmodel$global$calls)) {
    tcamglobe <- NULL
  } else{
    tcamglobe <- qpmodel$global$calls$r
  }


  #####  write components to file #####

  more <- function(x = '\n', sep = '\n', append = TRUE){
    cat(x, file = filename, sep = sep, append = append)
  }

  # initialize the file
  more('', sep = '', append = FALSE)

  # write GLOBAL
  if(length(tcamglobe)){
    more('$GLOBAL')
    more(tcamglobe)
    more()
  }

  # write parameters
  if(length(finalparams)){
    more('$PARAM')
    more(finalparams, sep = ', ')
    more()
  }

  # write compartments
  if(length(finalcmts)){
    more('$CMT ', sep = ' ')
    more(finalcmts, sep = ' ')
    more()
  }

  # write MAIN
  if(length(finalmain)){
    more('$MAIN')
    if(length(finalmain)) more(finalmain)
    if(length(tcammain)) more(tcammain)
    more('')
  }


  # write OMEGA
  if(length(omega)){
    more('$OMEGA @labels', sep = ' ')
    more(' ',sep = ' ')
    more(names(qpmodel$omega$calls), sep = ' ')
    more()
    more(omega,sep = ' ')
    more()
  }

  # write ODEs
  if(length(tcamode) | length(algebra) | length(ode)){
    more('$ODE')
    if(length(tcamode)) more(tcamode)
    if(length(algebra)) more(algebra)
    if(length(ode)) more(ode)
    more('')
  }

  # write SIGMA
  if(length(sigma)){
    more('$SIGMA @labels ', sep = ' ')
    more(' ', sep = ' ')
    more(names(qpmodel$sigma$calls),sep = ' ')
    more()
    more(sigma, sep = ' ')
    more()
  }

  # write TABLE
  if(length(table)){
    more('$TABLE')
    more(table)
    more('')
  }

  # write CAPTURE
  if(length(capture)){
    more('$CAPTURE ', sep = ' ')
    more(capture, sep = ' ')
    more('')
  }
  close(filemodel)
}

#' Coerce to NONMEM
#'
#' Coerces model to NONMEM format
#'
#' @param qpmodel class 'qpmodel'
#' @param filename output filename
#' @export
#'
modelCreate_NM <- function(
  qpmodel = qpModel(),
  filename = 'run1.mod'
){
  more <- function(x, sep = '\n')cat(x, file = filename, sep = sep, append = TRUE)

  ## Create model text file
  filemodelNM <- file(filename)

  ## start writing model code for NM file
  writeLines(c('$PROBLEM'), filemodelNM)
  more(';; 1. Based on:')
  more(';; 2. Description:')
  more(';; 3. Label:')
  more(';; 4. Structural model:')
  more(';; 5. Covariate model:')
  more(';; 6. Inter-individual variability:')
  more(';; 7. Inter-occasion variability:')
  more(';; 8. Residual variability:')
  more(';; 9. Estimation:')
  more('$INPUT')
  more('$DATA')
  more('$SUBROUTINE ADVAN13 TOL=12')

  ## create compartment list
  finalcmts <- names(qpmodel$ode$calls)
  NMcmts <- paste('COMP', '(', names(qpmodel$ode$calls), ')')

  ## create parameters using MU modeling
  theta_param <- paste(
    names(qpmodel$theta$calls),
    '=',
    'THETA(',
    c(1:length(qpmodel$theta$calls)),
    ')',
    sep = ''
  )
  mu_params <- paste(
    'MU_',
    c(1:length(qpmodel$theta$calls)),
    ' = ',
    names(qpmodel$theta$calls),
    sep = ''
  )
  final_params <- paste(names(qpmodel$param$calls), ' = ', qpmodel$param$calls)
  for (i in length(qpmodel$omega$calls):1) {
    final_params <- str_replace(
      final_params,
      names(qpmodel$omega$calls[i]),
      paste('ETA(', i, ')', sep = '')
    )
  }

  ## create ODEs
  states <- paste(
    'A(',
    1:length(names(qpmodel$ode$calls)),
    ')',
    ' = ',
    names(qpmodel$ode$calls),
    sep = ''
  )
  algebra <- paste(
    names(qpmodel$algebraic$calls),
    ' = ',
    qpmodel$algebraic$calls,
    sep = ''
  )
  des <- paste(
    'DADT(',
    1:length(names(qpmodel$ode$calls)),
    ')',
    ' = ',
    qpmodel$ode$calls,
    sep = ''
  )

  ### TCAM components

  ## MAIN & ODE
  if (is.null(qpmodel$custom$calls)) {
    tcamodenm <- NULL
  } else{
    tcamodenm <- qpmodel$custom$calls$nm.ode
  }

  ## GLOBAL
  if (is.null(qpmodel$global$calls)) {
    tcamglobenm <- NULL
  } else{
    tcamglobenm <- qpmodel$global$calls$nm
  }

  ## create ERROR block


  ## write CMTs section
  more('$MODEL',sep = '     ')
  more(NMcmts, sep = '  ')

  ## write PK block
  more()
  more('$PK')
  more(theta_param)
  more(mu_params)
  more(final_params)
  more()
  more(tcamglobenm)

  ## write DES block
  more()
  more('$DES')
  more(tcamodenm)
  more()
  more(states)
  more(algebra)
  more(des)
  more()

  ## write ERROR block
  more('$ERROR')

  more('Y = IPRED*(1+EPS(1)) + EPS(2)')
  more('W = SQRT(') # was just 'cat' in reference implementation
  more('IRES = DV - IPRED')
  more('IWRES = IRES/W')

  ## write THETAs
  more('$THETA')
  #more(qpmodel$theta$calls)
  more(unlist(qpmodel$theta$calls))

  ## write OMEGA
  more('$OMEGA')
  #more(qpmodel$omega$calls)
  more(unlist(qpmodel$omega$calls))

  ## write SIGMA
  more('$SIGMA')
  more(unlist(qpmodel$sigma$calls))

  ## write ESTIMATION
  more('$ESTIMATION METHOD=FOCE INTER NSIG=3 SIGL=9 MAXEVALS=9999 PRINT=10 NOABORT')

  ## write COVARIANCE
  more('$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S')

  ## write TABLE
  more('$TABLE')

  close(filemodelNM)

}

##### absorption models #####

#' Zero Order Absorption
#'
#' zero order absorption.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.abs.zero.order <- function(
  theta = eqns(tvk0 = 1),
  omega = eqns(etak0 = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0),
    algebraic = eqns(abs0 = k0),
    param = eqns(k0 = exp(tvk0 + etak0)),
    theta = theta,
    omega = omega,
    output = eqns(
      k0 = k0,
      etak0 = etak0
    )
  )
}

#' First Order Absorption
#'
#' First order absorption.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.abs.first.order <- function(
  theta = eqns(tvka = 1),
  omega = eqns(etaka = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0),
    algebraic = eqns(abs0 = ka * Aabs),
    param = eqns(ka = exp(tvka + etaka)),
    theta = theta,
    omega = omega,
    output = eqns(ka = ka,
                  etaka = etaka)
  )
}

#' Zero and First Order Absorption
#'
#' zero and first order absorption.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.abs.zero.and.first.order <- function(
  theta = eqns(tvka = 1, tvk0 = 1),
  omega = eqns(etaka = 0, etak0 = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0),
    algebraic = eqns(abs0 = +ka * Aabs +
                       k0),
    param = eqns(ka = exp(tvka + etaka),
                  k0 = exp(tvk0 + etak0)),
    theta = theta,
    omega = omega,
    output = eqns(
      ka = ka,
      k0 = k0,
      etaka = etaka,
      etak0 = etak0
    )
  )
}

#' Michaelis-Menten Absorption
#'
#' Michaelis-Menten absorption.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.abs.MM <- function(
  theta = eqns(tvamax = 1,tvka50 = 1),
  omega = eqns(etaka50 = 0, etaamax = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0),
    algebraic = eqns(abs0 = (amax * Aabs) / (ka50 + Aabs)),
    param = eqns(
      amax = exp(tvamax + etaamax),
      ka50 = exp(tvka50 + etaka50)
    ),
    theta = theta,
    omega = omega,
    output = eqns(
      amax = amax,
      ka50 = ka50,
      etaamax = etaamax,
      etaka50 = etaka50
    )
  )
}


#' Transit Compartment Absorption Model
#'
#' Transit compartment absorption model (TCAM).
#' @param algebraics algebraic equations
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.abs.TCAM <- function(
  algebraics = eqns(abs0 = inpt),
  theta  = eqns(tvmtt = 1, tvntr = 1, tvf1 = 0.1),
  omega = eqns(etamtt = 0, etantr = 0, etaf1 = 0)
){
  qpModel(
    param = eqns(
      MTT = exp(tvmtt + etamtt),
      NTR = exp(tvntr + etantr),
      LF1 = exp(tvf1 + etaf1)
    ),
    theta = theta,
    omega = omega,
    output = eqns(
      MTT = MTT,
      NTR = NTR,
      LF1 = LF1,
      etamtt = etamtt,
      etantr = etantr,
      etaf1 = etaf1
    ),
    global = eqns(
      r = '// initialize dose amt and dose time arrays
#define LEND 1000
#define MAXD 10
double damt[1000];
double dtime[1000];
int ndose=-1;

int counter;',
      nm = 'MAXD = 10
;; Initiate TCAM dose history
"IF(NEWIND.NE.2.OR.EVID.EQ.4) THEN
"  I = 1
"  DO WHILE (I.LE.(2*MAXD))
"    COM(I) = 0
"    I = I + 1
"  END DO
"END IF

;; Log a dose
"IF(AMT.NE.0.AND.(EVID.EQ.4.OR.EVID.EQ.1)) THEN
"  IF(COM(2*MAXD).EQ.0) THEN
"! We have not yet reached the maximum number of doses in history
"    I = 1
"    DO WHILE (I.LE.MAXD)
"      IF(COM(MAXD + I).EQ.0) THEN
"        COM(I)        = TIME
"        COM(MAXD + I) = AMT
"!        IF(ID.EQ.2051) WRITE(*,*) I, MAXD+I, COM(1), COM(2), COM(3), COM(4)
"	 I = MAXD+1
"      END IF
"      I = I + 1
"    END DO
"  ELSE
"! Maximum number of doses in history reached - need to reinit
"    I = 1
"    DO WHILE (I.LE.(MAXD-1))
"      COM(I)        = COM(I + 1)
"      COM(MAXD + I) = COM(MAXD + I + 1)
"    END DO
"    COM(MAXD)   = TIME
"    COM(2*MAXD) = AMT
"  END IF
"END IF

; TRANSIT compartment modeling
; Stirling approximation to gamma
X = 0.00001
;LOGF = 0.5*LOG(NTR*2*3.1415926+X) + NTR*LOG(NTR+X) - NTR + LOG(1+1/12/NTR+1/288/NTR/NTR+X)
LOGF =  GAMLN(NTR+1)

F1 = 0
X = 0.00001'
    ),
    custom = eqns(
      r.main = '// TCAM variables; be sure to define NTR and MTT
double KTR = (NTR+1)/MTT;
double DOSF1 = (LF1)/(1+LF1);
F_Acentral= 0; // set F1 to zero to avoid dosing into compartment.  were using TCAM here.
// Capture dose information
  if (NEWIND <= 1) {
    ndose = 0;
  }


  if(EVID == 1|EVID == 4){
    dtime[ndose] = TIME;
    damt[ndose] = self.amt*DOSF1;
    ndose++;
  }
                                    ',
      r.ode = 'double LOGF = lgamma(NTR+1);

double inpt = 0;
counter = ndose-1;

  // count down, not up, and only use last MAXD items
  while(counter>ndose-MAXD & counter>=0){
    double WTIME = SOLVERTIME - dtime[counter];
    double DOSI = damt[counter];
    double inp=0;
    if(WTIME>0.0 && DOSI>0.0){
      inp = exp(log(DOSI) + NTR*log(KTR*WTIME) + log(KTR) - KTR*WTIME - LOGF);
    }
    inpt += inp;
    counter--;

  } // end of while loop',
      nm.ode = 'INP = 0
 I = 1
 ;; Loop over all stored doses and sum INP
 "DO WHILE (I.LE.MAXD)
 "  IF(COM(MAXD + I).GT.0) THEN
 "    TSTAR = T - COM(I)
 "    DOSTC = COM(MAXD + I)
 "    IF(TSTAR.LT.0.001) TSTAR = 0.001
      INP = INP + EXP(LOG(DOSTC*DOSF1+X) + NTR*LOG(KTR*TSTAR+X) + LOG(KTR+X)- KTR*TSTAR - LOGF)
 "  END IF
 "  I = I + 1
 "END DO'
    ),
    algebraic = algebraics
  )



}

#####  elimination models #####

#' Elimination Model
#'
#' Elimination model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.elim <- function(
  theta = eqns(tvke = 1),
  omega = eqns(etake = 0)
){
  qpModel(
    algebraic = eqns(elim = -ke * Acentral),
    param = eqns(ke = exp(tvke + etake)),
    theta = eqns(tvke = 1),
    omega = eqns(etake = 0),
    output = eqns(ke = ke,
                  etake = etake)
  )
}

#' Clearance Model
#'
#' Clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.elim.clearance <- function(
  theta = eqns(tvCL = 1, tvV = 1),
  omega = eqns(etaCL = 0, etaV = 0)
){
  qpModel(
    algebraic = eqns(elim = -(CL / V) * Acentral),
    param = eqns(CL = exp(tvCL + etaCL),
                  V = exp(tvV + etaV)),
    theta = theta,
    omega = omega,
    output = eqns(
      CL = CL,
      V = V,
      etaCL = etaCL,
      etaV = etaV
    )
  )
}

#' Michaelis-Menten Clearance
#'
#' Micahaelis-Menten clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @export
#' @family models
qp.elim.MM <- function(
  theta = eqns(tvkmax = 1, tvkm50 = 1),
  omega = eqns(etakmax = 0, etakm50 = 0)
){
  qpModel(
    algebraic = eqns(elim = -(kmax * Acentral) / (km50 + Acentral)),
    param = eqns(
      kmax = exp(tvkmax + etakmax),
      km50 = exp(tvkm50 + etakm50)
    ),
    theta = theta,
    omega = omega,
    output = eqns(
      kmax = kmax,
      km50 = km50,
      etakmax = etakmax,
      etakm50 = etakm50
    )
  )
}

#####  compartmental models #####
#' One Compartment IV Model (with Rate)
#'
#' One compartment IV model estimating infusion rate.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.onecpt.iv.rate <- function(
  theta = eqns(tvke = 1, tvV = 1),
  omega = eqns(etake = 0, etaV = 0),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(Acentral = abs0 - elim),
    algebraic = eqns(abs0 = 0,
                     elim = ke * Acentral),
    param = eqns(ke = exp(tvke + etake),
                  V = exp(tvV + etaV)),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      ke = ke,
      V = V,
      etake = etake,
      etaV = etaV
    ),
    custom = eqns(),
    global = eqns()
  )
}

#' One Compartment IV Clearance Model
#'
#' One compartment IV clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.onecpt.iv.cl <- function(
  theta = eqns(tvcl = 1, tvV = 1),
  omega = eqns(etacl = 0, etaV = 0),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(Acentral = abs0 - elim),
    algebraic = eqns(abs0 = 0,
                     elim = (cl / V) * Acentral),
    param = eqns(cl = exp(tvcl + etake),
                  V = exp(tvV + etaV)),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      cl = cl,
      V = V,
      etacl = etacl,
      etaV = etaV
    ),
    custom = eqns(),
    global = eqns()
  )
}

#' One Compartment Oral Fixed Rate Clearance
#'
#' One compartment oral fixed rate clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.onecpt.oral.rate <- function(
  theta = eqns(tvke = 1,tvV = 1, tvka = 1),
  omega = eqns(etake = 0, etaV = 0, tvka = 0),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0,
                Acentral = abs0 - elim),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = ke * Acentral),
    param = eqns(
      ke = exp(tvke + etake),
      ka = exp(tvka + etaka),
      V = exp(tvV + etaV)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      ke = ke,
      V = V,
      ka = ka,
      etake = etake,
      etaV = etaV
    )
  )
}

#' One Compartment Oral Clearance Model
#'
#' One compartment oral clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.onecpt.oral.cl <- function(
  theta = eqns(tvcl = 1, tvV = 1, tvka = 1),
  omega = eqns(etacl = 0, etaV = 0, etaka = 0),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(Aabs = -abs0,
                Acentral = abs0 - elim),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = (cl / V) * Acentral),
    param = eqns(
      cl = exp(tvcl + etacl),
      ka = exp(tvka + etaka),
      V = exp(tvV + etaV)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      cl = cl,
      V = V,
      ka = ka,
      etacl = etacl,
      etaV = etaV,
      etaka = etaka
    )
  )
}

#' Two Compartment IV Rate Model
#'
#' Two compartment IV model estimating infusion rate.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.twocpt.iv.rate <- function(
  theta = eqns(
    tvke = 1,
    tvV = 1,
    tvk12 = 1,
    tvk21 = 1
  ),
  omega = eqns(
    etake = 0,
    etaV = 0,
    etak12 = 0,
    etak21 = 0
  ),
  sigma = eqns(
    add = 0, prop = 0)
){
    qpModel(
      ode = eqns(
        Acentral = abs0 - elim - k12 * Acentral + k21 * Aperiph,
        Aperiph = k12 * Acentral - k21 * Aperiph
      ),
      algebraic = eqns(abs0 = 0,
                       elim = ke * Acentral),
      param = eqns(
        ke = exp(tvke + etake),
        V = exp(tvV + etaV),
        k12 = exp(tvk12 + etak12),
        k21 = exp(tvk21 + etak21)
      ),
      theta = theta,
      omega = omega,
      sigma = sigma,
      observe = eqns(IPRED = Acentral / V,
                     DV = IPRED * (1 + prop) + add),
      output = eqns(
        IPRED = IPRED,
        DV = DV,
        ke = ke,
        V = V,
        k12 = k12,
        k21 = k21,
        etake = etake,
        etaV = etaV,
        etak12 = etak12,
        etak21 = etak21
      )
    )
}

#' Two Compartment IV Clearance Model
#'
#' Two compartment IV clearance model
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.twocpt.iv.cl <- function(
  theta = eqns(
    tvCL = 1,
    tvV = 1,
    tvQ = 1,
    tvV2 = 1
  ),
  omega = eqns(
    etaCL = 0,
    etaV = 0,
    etaQ = 0,
    etaV2 = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Acentral = abs0 - elim - Q * (Acentral / V - Aperiph / V2),
      Aperiph = Q * (Acentral / V - Aperiph / V2)
    ),
    algebraic = eqns(abs0 = 0,
                     elim = (CL / V) * Acentral),
    param = eqns(
      CL = exp(tvCL + etaCL),
      V = exp(tvV + etaV),
      Q = exp(tvQ + etaQ),
      V2 = exp(tvV2 + etaV2)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      CL = CL,
      V = V,
      Q = Q,
      V2 = V2,
      etaCL = etaCL,
      etaV = etaV,
      etaQ = etaQ,
      etaV2 = etaV2
    )
  )
}

#' Two Compartment Oral Rate Model
#'
#' Two compartment oral rate model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.twocpt.oral.rate <- function(
  theta = eqns(
    tvke = 1,
    tvV = 1,
    tvk12 = 1,
    tvk21 = 1,
    tvka = 1
  ),
  omega = eqns(
    etake = 0,
    etaV = 0,
    etak12 = 0,
    etak21 = 0,
    etaka = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Aabs = -abs0,
      Acentral = abs0 - elim - k12 * Acentral + k21 *
        Aperiph,
      Aperiph = k12 * Acentral - k21 *
        Aperiph
    ),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = ke * Acentral),
    param = eqns(
      ke = exp(tvke + etake),
      V = exp(tvV + etaV),
      k12 = exp(tvk12 + etak12),
      k21 = exp(tvk21 + etak21),
      ka = exp(tvka + etaka)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      ke = ke,
      V = V,
      k12 = k12,
      k21 = k21,
      ka = ka,
      etake = etake,
      etaV = etaV,
      etak12 = etak12,
      etak21 = etak21,
      etaka = etaka
    )
  )
}

#' Two Compartment Oral Clearance Model
#'
#' Two compartment oral clearance model
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.twocpt.oral.cl <- function(
  theta = eqns(
    tvCL = 1,
    tvV = 1,
    tvQ = 1,
    tvV2 = 1,
    tvka = 1
  ),
  omega = eqns(
    etaCL = 0,
    etaV = 0,
    etaQ = 0,
    etaV2 = 0,
    etaka = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Aabs = -abs0,
      Acentral = abs0 - elim - Q * (Acentral / V - Aperiph /
                                      V2),
      Aperiph = Q * (Acentral / V - Aperiph / V2)
    ),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = (CL / V) * Acentral),
    param = eqns(
      CL = exp(tvCL + etaCL),
      V = exp(tvV + etaV),
      Q = exp(tvQ + etaQ),
      V2 = exp(tvV2 + etaV2),
      ka = exp(tvka + etaka)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      CL = CL,
      V = V,
      Q = Q,
      V2 = V2,
      ka = ka,
      etaCL = etaCL,
      etaV = etaV,
      etaQ = etaQ,
      etaV2 = etaV2,
      etaka = etaka
    )
  )
}

#' Three Compartment IV Rate Model
#'
#' Three compartment IV rate model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.threecpt.iv.rate <- function(
  theta = eqns(
    tvke = 1,
    tvV = 1,
    tvk12 = 1,
    tvk21 = 1,
    tvk13 = 1,
    tvk31 = 1
  ),
  omega = eqns(
    etake = 0,
    etaV = 0,
    etak12 = 0,
    etak21 = 0,
    etak13 = 0,
    etak31 = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Acentral = abs0 - elim - k12 * Acentral + k21 * Aperiph - k13 * Acentral + k31 *
        A2periph,
      Aperiph = k12 * Acentral - k21 * Aperiph,
      A2periph = k13 * Acentral - k31 * A2periph
    ),
    algebraic = eqns(abs0 = 0,
                     elim = ke * Acentral),
    param = eqns(
      ke = exp(tvke + etake),
      V = exp(tvV + etaV),
      k12 = exp(tvk12 + etak12),
      k21 = exp(tvk21 + etak21),
      k13 = exp(tvk13 + etak13),
      k31 = exp(tvk31 + etak31)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      ke = ke,
      V = V,
      k12 = k12,
      k21 = k21,
      k13 = k13,
      k31 = k31,
      etake = etake,
      etaV = etaV,
      etak12 = etak12,
      etak21 = etak21
    )
  )
}

#' Three Compartment IV Clearance Model
#'
#' Three compartment IV clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.threecpt.iv.cl <- function(
  theta = eqns(
    tvCL = 1,
    tvV = 1,
    tvQ = 1,
    tvV2 = 1,
    tvQ2 = 1,
    tvV3 = 1
  ),
  omega = eqns(
    etaCL = 0,
    etaV = 0,
    etaQ = 0,
    etaV2 = 0,
    etaQ2 = 0,
    etaV3 = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Acentral = abs0 - elim - Q * (Acentral / V - Aperiph / V2) - Q2 * (Acentral /
                                                                           V - A2periph / V3),
      Aperiph = Q * (Acentral / V - Aperiph / V2),
      A2periph = Q2 * (Acentral / V - A2periph / V3)
    ),
    algebraic = eqns(abs0 = 0,
                     elim = (CL / V) * Acentral),
    param = eqns(
      CL = exp(tvCL + etaCL),
      V = exp(tvV + etaV),
      Q = exp(tvQ + etaQ),
      V2 = exp(tvV2 + etaV2),
      Q2 = exp(tvQ2 + etaQ2),
      V3 = exp(tvV3 + etaV3)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      CL = CL,
      V = V,
      Q = Q,
      V2 = V2,
      Q2 = Q2,
      V3 = V3,
      etaCL = etaCL,
      etaV = etaV,
      etaQ = etaQ,
      etaV2 = etaV2
    )
  )
}

#' Three compartment Oral Rate Model
#'
#' Three compartment oral rate model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.threecpt.oral.rate <- function(
  theta = eqns(
    tvke = 1,
    tvV = 1,
    tvk12 = 1,
    tvk21 = 1,
    tvk13 = 1,
    tvk31 = 1,
    tvka = 1
  ),
  omega = eqns(
    etake = 0,
    etaV = 0,
    etak12 = 0,
    etak21 = 0,
    etak13 = 0,
    etak31 = 0,
    etaka = 0
  ),
  sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Aabs = -abs0,
      Acentral = abs0 - elim - k12 * Acentral + k21 * Aperiph - k13 *
        Acentral + k31 * A2periph,
      Aperiph = k12 * Acentral - k21 * Aperiph,
      A2periph = k13 * Acentral - k31 * A2periph
    ),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = ke * Acentral),
    param = eqns(
      ke = exp(tvke + etake),
      V = exp(tvV + etaV),
      k12 = exp(tvk12 + etak12),
      k21 = exp(tvk21 + etak21),
      k13 = exp(tvk13 + etak13),
      k31 = exp(tvk31 + etak31),
      ka = exp(tvka + etaka)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      ke = ke,
      V = V,
      k12 = k12,
      k21 = k21,
      k13 = k13,
      k31 = k31,
      ka = ka,
      etake = etake,
      etaV = etaV,
      etak12 = etak12,
      etak21 = etak21,
      etak13 = etak13,
      etak31 = etak31,
      etaka = etaka
    )
  )
}

#' Three Compartment Oral Clearance Model
#'
#' Three compartment oral clearance model.
#' @param theta fixed effects
#' @param omega second level random effects
#' @param sigma first level random effects
#' @export
#' @family models
qp.threecpt.oral.cl <- function(
  theta = eqns(
  tvCL = 1,
  tvV = 1,
  tvQ = 1,
  tvV2 = 1,
  tvQ2 = 1,
  tvV3 = 1,
  tvka = 1
),
omega = eqns(
  etaCL = 0,
  etaV = 0,
  etaQ = 0,
  etaV2 = 0,
  etaQ2 = 0,
  etaV3 = 0,
  etaka = 0
),
sigma = eqns(add = 0, prop = 0)
){
  qpModel(
    ode = eqns(
      Aabs = -abs0,
      Acentral = abs0 - elim - Q * (Acentral / V - Aperiph /
                                      V2) - Q2 * (Acentral / V - A2periph / V3),
      Aperiph = Q * (Acentral / V - Aperiph / V2),
      A2periph = Q2 * (Acentral / V - A2periph / V3)
    ),
    algebraic = eqns(abs0 = ka * Aabs,
                     elim = (CL / V) * Acentral),
    param = eqns(
      CL = exp(tvCL + etaCL),
      V = exp(tvV + etaV),
      Q = exp(tvQ + etaQ),
      V2 = exp(tvV2 + etaV2),
      Q2 = exp(tvQ2 + etaQ2),
      V3 = exp(tvV3 + etaV3),
      ka = exp(tvka + etaka)
    ),
    theta = theta,
    omega = omega,
    sigma = sigma,
    observe = eqns(IPRED = Acentral / V,
                   DV = IPRED * (1 + prop) + add),
    output = eqns(
      IPRED = IPRED,
      DV = DV,
      CL = CL,
      V = V,
      Q = Q,
      V2 = V2,
      Q2 = Q2,
      V3 = V3,
      ka = ka,
      etaCL = etaCL,
      etaV = etaV,
      etaQ = etaQ,
      etaV2 = etaV2,
      etaka = etaka
    )
  )
}


##########   update functions ##########

#' Update Thetas
#'
#' Updates theta.
#' @export
#' @family updates
#' @param theta fixed effects
theta <- function(theta = eqns()) {
  qpModel(theta = theta)
}

#' Update Omegas
#'
#' Updates omega.
#' @export
#' @family updates
#' @param omega second level random effects
omega <- function(omega = eqns()) {
  qpModel(omega = omega)
}

#' Update Sigmas
#'
#' Updates sigma.
#' @param sigma first level random effects
#' @export
#' @family updates
sigma <- function(sigma = eqns()) {
  qpModel(sigma = sigma)
}

#' Update Parameter Definitions
#'
#' Updates parameter definitions.
#' @export
#' @family updates
#' @param param parameters
parameter <- function(param = eqns()) {
  qpModel(param = param)
}

#' Update ODE
#'
#' Updates ODEs.
#' @export
#' @family updates
#' @param ode ondinary differential equations
ode <- function(ode = eqns()) {
  qpModel(ode = ode)
}
