#' Create a qpModel

#' Creates an object of class 'qpmodel', convertible
#' to NONMEM code or mrgsolve code.
#'
#' @param odes ordinary differential equations
#' @param algebraic algebraic equations
#' @param params parameters
#' @param omega second level random effects
#' @param sigma first level random effects
#' @param theta fixed effects
#' @param observe items to include in output
#' @param output other output
#' @param global global equations
#' @param custom custom equations
#' @export
#' @examples
#' qpModel()
#'
qpModel = function(odes = eqns(),
                   algebraic = eqns(),
                   params = eqns(),
                   omega = eqns(),
                   sigma = eqns(),
                   theta = eqns(),
                   observe = eqns(),
                   output = eqns(),
                   global = eqns(),
                   custom = eqns()
                   ){


  me = list(
    odes = odes,
    algebraic = algebraic,
    params = params,
    omega = omega,
    sigma = sigma,
    theta = theta,
    observe = observe,
    output = output,
    global = global,
    custom = custom
  )

  class(me) = "qpmodel"
  return(me)

}

eqns = function(...){
  exprns = enexprs(...)
  structure(list(calls = exprns),class="eqns")
}

eqns.make = function(...){
  structure(list(calls=arg),class="eqns")
}

eqn.make = function(arg){
  # arg is a class::call

  # this creates a string from the call for operations
  instring=deparse(arg)
  structure(list(call=arg,
                 instring=instring,
                 vars=all.vars(arg)),
            class="eqn")
}

eqn = function(arg){
  # arg is an expression
  exprn = enexpr(arg)

  eqn.make(exprn)
  #instring = deparse(arg)
  #structure(list(call = exprn,instring = instring),class="eqn")
}

eqn_ = function(arg){
  # arg is a expression

  # this converts arg into a call
  arg$new = arg
  arg = arg[[1]]
  eqn.make(arg)
  #structure(list(call = arg),class="eqn")
}

eqns_ = function(...){
  # arg is a expression
  exprns = enexprs(...)

  # this converts arg into a call
  exprns$new = exprns
  exprns = exprns[length(match.call)-1]
  #eqn.make(exprns)
  structure(list(calls = exprns),class="eqns")
}

`+.eqn` <- function(eq1,eq2){

  neweq = parse(text=paste("(",eq1$instring,")+(",eq2$instring,")"))
  eqn_(neweq)

}

`-.eqn` <- function(eq1,eq2){

  neweq = parse(text=paste("(",eq1$instring,")-(",eq2$instring,")"))
  eqn_(neweq)

}

`*.eqn` <- function(eq1,eq2){

  neweq = parse(text=paste("(",eq1$instring,")*(",eq2$instring,")"))
  eqn_(neweq)

}

`/.eqn` <- function(eq1,eq2){

  neweq = parse(text=paste("(",eq1$instring,")/(",eq2$instring,")"))
  eqn_(neweq)


}

`^.eqn` <- function(eq1,eq2){

  neweq = parse(text=paste("(",eq1$instring,")^(",eq2$instring,")"))
  eqn_(neweq)

}

`+.eqns` <- function(eqns1,eqns2){

  # set up new eqns object
  neweqns = eqns()

  # find calls in both eqn structures
  interx = intersect(names(eqns1$calls),names(eqns2$calls))

  # find location of duplicates in eqns1
  locate = names(eqns1$calls) %in% interx

  # create addition of two objects
  neweqns$calls = c(eqns1$calls[!locate],eqns2$calls)

  structure(list(calls=neweqns$calls),class="eqns")

}

`*.eqns` <- function(eqns1,eqns2){

  # set up new eqns object
  neweqns = eqns()

  # create addition of two objects
  neweqns$calls = c(eqns1$calls,eqns2$calls)

  structure(list(calls=neweqns$calls),class="eqns")

}

`+.qpmodel` = function(qp1,qp2){

  #odesnew = append(qp1$odes$calls,qp2$odes$calls)
  #algebraicnew = append(qp1$algebraic$calls,qp2$algebraic$calls)
  #paramsnew = append(qp1$params$calls,qp2$params$calls)
  #thetanew = append(qp1$theta$calls,qp2$theta$calls)
  #omeganew = append(qp1$omega$calls,qp2$omega$calls)
  #sigmanew = append(qp1$sigma$calls,qp2$sigma$calls)
  qpnew = qpModel()

  qpnew$odes = qp1$odes+qp2$odes
  qpnew$algebraic = qp1$algebraic+qp2$algebraic
  qpnew$params = qp1$params+qp2$params
  qpnew$theta = qp1$theta+qp2$theta
  qpnew$omega = qp1$omega+qp2$omega
  qpnew$sigma = qp1$sigma+qp2$sigma
  qpnew$observe = qp1$observe + qp2$observe
  qpnew$output = qp1$output + qp2$output
  qpnew$global = qp1$global + qp2$global
  qpnew$custom = qp1$custom + qp2$custom

  #odesnew = eqn_(odesnew)

  #qp1$odes = qp1$odes + qp2$odes
  #qp1$algebraic = qp1$algebraic + qp2$algebraic
  #qp1$params = qp1$params + qp2$params
  #qp1$theta = qp1$theta + qp2$theta
  #qp1$omega = qp1$omega + qp2$omega
  #qp1$sigma = qp1$sigma + qp2$omega

  qpModel(qpnew$odes,
          qpnew$algebraic,
          qpnew$params,
          qpnew$omega,
          qpnew$sigma,
          qpnew$theta,
          qpnew$observe,
          qpnew$output,
          qpnew$global,
          qpnew$custom)
}

`*.qpmodel` = function(qp1,qp2){

  #odesnew = append(qp1$odes$calls,qp2$odes$calls)
  #algebraicnew = append(qp1$algebraic$calls,qp2$algebraic$calls)
  #paramsnew = append(qp1$params$calls,qp2$params$calls)
  #thetanew = append(qp1$theta$calls,qp2$theta$calls)
  #omeganew = append(qp1$omega$calls,qp2$omega$calls)
  #sigmanew = append(qp1$sigma$calls,qp2$sigma$calls)
  qpnew = qpModel()

  qpnew$odes = qp1$odes*qp2$odes
  qpnew$algebraic = qp1$algebraic*qp2$algebraic
  qpnew$params = qp1$params*qp2$params
  qpnew$theta = qp1$theta*qp2$theta
  qpnew$omega = qp1$omega*qp2$omega
  qpnew$sigma = qp1$sigma*qp2$sigma
  qpnew$observe = qp1$observe * qp2$observe
  qpnew$output = qp1$output * qp2$output
  qpnew$global = qp1$global * qp2$global
  qpnew$custom = qp1$custom * qp2$custom

  #odesnew = eqn_(odesnew)

  #qp1$odes = qp1$odes + qp2$odes
  #qp1$algebraic = qp1$algebraic + qp2$algebraic
  #qp1$params = qp1$params + qp2$params
  #qp1$theta = qp1$theta + qp2$theta
  #qp1$omega = qp1$omega + qp2$omega
  #qp1$sigma = qp1$sigma + qp2$omega

  qpModel(qpnew$odes,
          qpnew$algebraic,
          qpnew$params,
          qpnew$omega,
          qpnew$sigma,
          qpnew$theta,
          qpnew$observe,
          qpnew$output,
          qpnew$global,
          qpnew$custom)
}

modelView = function(qpmodel = qpModel()){

  # write out odes
  #print("ODEs:")
  #print(paste(names(qpmodel$odes$calls)," = ",qpmodel$odes$calls))
  odes = data.frame("ODEs" = paste(names(qpmodel$odes$calls)," = ",qpmodel$odes$calls))

  # write out relationships
  #print("Algebraic eqns:")
  #print(paste(names(qpmodel$algebraic$calls)," = ",qpmodel$algebraic$calls ))
  algebra = data.frame("Algebraic" = paste(names(qpmodel$algebraic$calls)," = ",qpmodel$algebraic$calls))

  # write out parameters
  #print("Parameters:")
  #print(paste(names(qpmodel$params$calls)," = ",qpmodel$params$calls ))
  params = data.frame("Parameters" = paste(names(qpmodel$params$calls)," = ",qpmodel$params$calls))

  # write out fixed effect values
  #print("Fixed Effect values:")
  #print(paste(names(qpmodel$theta$calls)," = ",qpmodel$theta$calls ))
  theta = data.frame("Theta" = paste(names(qpmodel$theta$calls)," = ",qpmodel$theta$calls))

  # write out omega values
  #print("Omega:")
  #print(paste(names(qpmodel$omega$calls)," = ",qpmodel$omega$calls ))
  omega = data.frame("Omega" = paste(names(qpmodel$omega$calls)," = ",qpmodel$omega$calls))

  # write out observation info
  #print("Observations:")
  #print(paste(names(qpmodel$observe$calls)," = ",qpmodel$observe$calls ))
  observe = data.frame("Observe" = paste(names(qpmodel$observe$calls)," = ",qpmodel$observe$calls))

  # write out error values
  #print("Residual Error:")
  #print(paste(names(qpmodel$sigma$calls)," = ",qpmodel$sigma$calls ))
  sigma = data.frame("Sigma" = paste(names(qpmodel$sigma$calls)," = ",qpmodel$sigma$calls))

  # write out output variables
  output = data.frame("Output" = paste(qpmodel$output$calls))

  preview = list("ODEs" = odes,
                 "Algebraic" = algebra,
                 "Params" = params,
                 "Theta" = theta,
                 "Omega" = omega,
                 "Observe" = observe,
                 "Sigma" = sigma,
                 "Output" = output)

  return(preview)
}

modelCreate_R = function(qpmodel = qpModel(),filename = "model.txt"){

  ## Create model text file
  filemodel = file(filename)

  ## start writing model code for mrgsolve
  writeLines(c("modelfile =  "),filemodel)

  ## create parameter list
  finalparams = paste(names(qpmodel$theta$calls),"=",qpmodel$theta$calls)

  ## create compartment list
  finalcmts = names(qpmodel$odes$calls)

  ## create main
  finalmain = paste("double",names(qpmodel$params$calls),"=",qpmodel$params$calls,";")

  ## create ODEs
  odes = paste("dxdt_",names(qpmodel$odes$calls)," = ",qpmodel$odes$calls,";",sep="")

  ## create Algebraic eqns
  algebra = paste("double",names(qpmodel$algebraic$calls),"=",qpmodel$algebraic$calls,";")

  ## create SIGMA
  sigma = paste(qpmodel$sigma$calls)

  ## create OMEGA
  omega = paste(qpmodel$omega$calls)

  ## create observation equations
  if(is.null(qpmodel$observe$call)){
    table = NULL
  } else{
    table = paste("double",names(qpmodel$observe$calls)," = ",qpmodel$observe$calls,";")
  }

  ## create table output
  capture = paste(names(qpmodel$output$calls))

  ### TCAM components

  ## MAIN & ODE
  if(is.null(qpmodel$custom$calls)){
    tcammain = NULL
    tcamode = NULL
  } else{
    tcammain = qpmodel$custom$calls$r.main
    tcamode = qpmodel$custom$calls$r.ode
  }

  ## GLOBAL
  if(is.null(qpmodel$global$calls)){
    tcamglobe = NULL
  } else{
    tcamglobe = qpmodel$global$calls$r
  }


  #####  write components to file #####

  # write GLOBAL
  cat("$GLOBAL", file = filename, sep = "\n", append=TRUE)
  cat(tcamglobe, file = filename, sep = "\n", append=TRUE)

  # write parameters
  cat("", file = filename, sep = "\n", append=TRUE)
  cat("$PARAM", file = filename, sep = "\n", append=TRUE)
  cat(finalparams,file = filename, sep = ", ", append=TRUE)

  # write compartments
  cat("", file = filename, sep = "\n", append=TRUE)
  cat("$CMT ", file = filename, sep = " ", append=TRUE)
  cat(finalcmts,file = filename, sep=" ",append=TRUE)

  # write MAIN
  cat("", file = filename, sep="\n", append = TRUE)
  cat("$MAIN", file = filename, sep = "\n", append=TRUE)
  cat(finalmain, file = filename, sep = "\n", append = TRUE)
  cat(tcammain, file = filename, sep = "\n", append = TRUE)

  # write OMEGA
  cat("$OMEGA @labels", file = filename, sep = " ", append = TRUE)
  cat(" ", file = filename, sep = " ", append=TRUE)
  cat(names(qpmodel$omega$calls), file = filename, sep=" ", append=TRUE)
  cat("", file = filename, sep="\n", append=TRUE )
  cat(omega,file = filename,sep = " ", append=TRUE)

  # write ODEs
  cat("", file = filename, sep = "\n", append = TRUE)
  cat("$ODE", file = filename, sep = "\n", append = TRUE)
  cat(tcamode, file = filename, sep = "\n", append = TRUE)
  cat(algebra,file = filename, sep = "\n", append = TRUE)
  cat(odes, file = filename, sep = "\n", append = TRUE)

  # write SIGMA
  cat("$SIGMA @labels ", file = filename, sep = " ", append=TRUE)
  cat(" ", file = filename, sep = " ", append=TRUE)
  cat(names(qpmodel$sigma$calls),file = filename, sep = " ", append = TRUE)
  cat("", file = filename, sep="\n", append = TRUE)
  cat(sigma,file = filename, sep = " ", append=TRUE)

  # write TABLE
  cat("", file = filename, sep = "\n", append=TRUE)
  cat("$TABLE", file = filename, sep = "\n", append = TRUE)
  cat(table, file = filename, sep = "\n", append = TRUE)

  # write CAPTURE
  cat("$CAPTURE ", file = filename, sep=" ", append = TRUE)
  cat(capture,file = filename,sep = " ", append=TRUE)
  cat("",file = filename, sep = "\n", append=TRUE)


  close(filemodel)


}

modelCreate_NM = function(qpmodel = qpModel(),filename = "run1.mod"){

  ## Create model text file
  filemodelNM = file(filename)

  ## start writing model code for NM file
  writeLines(c("$PROBLEM"),filemodelNM)
  cat(";; 1. Based on:", file=filename,sep = "\n", append=TRUE)
  cat(";; 2. Description:", file=filename,sep = "\n", append=TRUE)
  cat(";; 3. Label:", file=filename,sep = "\n", append=TRUE)
  cat(";; 4. Structural model:", file=filename,sep = "\n", append=TRUE)
  cat(";; 5. Covariate model:", file=filename,sep = "\n", append=TRUE)
  cat(";; 6. Inter-individual variability:", file=filename,sep = "\n", append=TRUE)
  cat(";; 7. Inter-occasion variability:", file=filename,sep = "\n", append=TRUE)
  cat(";; 8. Residual variability:", file=filename,sep = "\n", append=TRUE)
  cat(";; 9. Estimation:", file=filename,sep = "\n", append=TRUE)
  cat("$INPUT", file=filename,sep = "\n", append=TRUE)
  cat("$DATA", file=filename,sep = "\n", append=TRUE)
  cat("$SUBROUTINE ADVAN13 TOL=12", file=filename,sep = "\n", append=TRUE)

  ## create compartment list
  finalcmts = names(qpmodel$odes$calls)
  NMcmts = paste("COMP","(",names(qpmodel$odes$calls),")")

  ## create parameters using MU modeling
  theta_param = paste(names(qpmodel$theta$calls),"=","THETA(",c(1:length(qpmodel$theta$calls)),")",sep="")
  mu_params = paste("MU_",c(1:length(qpmodel$theta$calls))," = ",names(qpmodel$theta$calls),sep="")
  final_params = paste(names(qpmodel$params$calls)," = ",qpmodel$params$calls)
  for(i in length(qpmodel$omega$calls):1){
    final_params = str_replace(final_params,names(qpmodel$omega$calls[i]),paste("ETA(",i,")",sep=""))
  }

  ## create ODEs
  states = paste("A(",1:length(names(qpmodel$odes$calls)),")"," = ",names(qpmodel$odes$calls),sep="")
  algebra = paste(names(qpmodel$algebraic$calls)," = ",qpmodel$algebraic$calls,sep="")
  des = paste("DADT(",1:length(names(qpmodel$odes$calls)),")"," = ",qpmodel$odes$calls,sep="")

  ### TCAM components

  ## MAIN & ODE
  if(is.null(qpmodel$custom$calls)){
    tcamodenm = NULL
  } else{
    tcamodenm = qpmodel$custom$calls$nm.ode
  }

  ## GLOBAL
  if(is.null(qpmodel$global$calls)){
    tcamglobenm = NULL
  } else{
    tcamglobenm = qpmodel$global$calls$nm
  }

  ## create ERROR block


  ## write CMTs section
  cat("$MODEL", file = filename, sep="     ", append=TRUE)
  cat(NMcmts, file = filename, sep = "  ", append=TRUE)

  ## write PK block
  cat("", file = filename, sep="\n", append = TRUE)
  cat("$PK", file = filename, sep = "\n", append=TRUE)
  cat(theta_param, file = filename, sep = "\n", append=TRUE)
  cat(mu_params, file = filename, sep = "\n", append=TRUE)
  cat(final_params, file = filename, sep = "\n", append=TRUE)
  cat("",file=filename, sep="\n", append=TRUE)
  cat(tcamglobenm, file = filename, sep = "\n", append=TRUE)



  ## write DES block
  cat("", file = filename, sep = "\n", append=TRUE)
  cat("$DES", file = filename, sep = "\n", append = TRUE)
  cat(tcamodenm,file = filename, sep = "\n", append = TRUE)
  cat("", file = filename, sep ="\n", append = TRUE)
  cat(states, file =filename, sep = "\n", append=TRUE)
  cat(algebra, file = filename, sep = "\n", append=TRUE)
  cat(des, file = filename, sep = "\n", append=TRUE)
  cat("", file = filename, sep = "\n", append=TRUE)

  ## write ERROR block
  cat("$ERROR", file = filename, sep = "\n", append=TRUE)

  cat("Y = IPRED*(1+EPS(1)) + EPS(2)", file = filename, sep = "\n", append=TRUE)
  cat("W = SQRT(")
  cat("IRES = DV - IPRED",file = filename, sep = "\n", append=TRUE)
  cat("IWRES = IRES/W")

  ## write THETAs
  cat("$THETA",file = filename, sep = "\n", append = TRUE)
  #cat(qpmodel$theta$calls, file = filename, sep = "\n", append=TRUE)
  cat(unlist(qpmodel$theta$calls),file = filename,sep = "\n", append=TRUE)

  ## write OMEGA
  cat("$OMEGA", file = filename, sep = "\n", append=TRUE)
  #cat(qpmodel$omega$calls, file  = filename, sep = "\n", append=TRUE)
  cat(unlist(qpmodel$omega$calls),file = filename,sep = "\n", append=TRUE)

  ## write SIGMA
  cat("$SIGMA", file = filename, sep = "\n", append=TRUE)
  cat(unlist(qpmodel$sigma$calls), file = filename, sep = "\n", append=TRUE)

  ## write ESTIMATION
  cat("$ESTIMATION METHOD=FOCE INTER NSIG=3 SIGL=9 MAXEVALS=9999 PRINT=10 NOABORT", file = filename, sep = "\n", append=TRUE)

  ## write COVARIANCE
  cat("$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S", file = filename, sep = "\n", append=TRUE)

  ## write TABLE
  cat("$TABLE", file = filename, sep = "\n", append=TRUE)

  close(filemodelNM)

}

##### absorption models #####

# zero order absorption
qp.abs.zero.order = function(thetas = eqns(tvk0 = 1),
                             omegas = eqns(etak0 = 0)){
  qpModel(odes = eqns(Aabs = -abs0),
          algebraic = eqns(abs0 = k0),
          params = eqns(k0 = exp(tvk0 + etak0)),
          theta = thetas,
          omega = omegas,
          output = eqns(k0=k0,
                        etak0=etak0)
          )
}
# first order absorption
qp.abs.first.order = function(thetas = eqns(tvka = 1),
                              omegas = eqns(etaka = 0)){
  qpModel(odes = eqns(Aabs = -abs0),
         algebraic = eqns(abs0 = ka*Aabs),
         params = eqns(ka = exp(tvka + etaka)),
         theta = thetas,
         omega = omegas,
         output = eqns(ka=ka,
                       etaka=etaka)
         )
}

# zero and first order absorption
qp.abs.zero.and.first.order = function(thetas = eqns(tvka = 1,
                                                     tvk0 = 1),
                                       omegas = eqns(etaka = 0,
                                                     etak0 = 0)){
  qpModel(odes = eqns(Aabs = -abs0),
                                      algebraic = eqns(abs0 = +ka*Aabs+k0),
                                      params = eqns(ka = exp(tvka + etaka),
                                                    k0 = exp(tvk0 + etak0)),
                                      theta = thetas,
                                      omega = omegas,
                                      output = eqns(ka=ka,
                                                    k0=k0,
                                                    etaka=etaka,
                                                    etak0=etak0)
                                      )
}

# michaelis menton absorption
qp.abs.MM = function(thetas = eqns(tvamax = 1,
                                   tvka50 = 1),
                     omegas = eqns(etaka50 = 0,
                                   etaamax = 0)){
  qpModel(odes = eqns(Aabs = -abs0),
          algebraic = eqns(abs0 = (amax*Aabs)/(ka50+Aabs)),
          params = eqns(amax = exp(tvamax + etaamax),
                        ka50 = exp(tvka50 + etaka50)),
          theta = thetas,
          omega = omegas,
          output = eqns(amax=amax,
                        ka50=ka50,
                        etaamax=etaamax,
                        etaka50=etaka50)
          )
}


# TCAM
qp.abs.TCAM = function(algebraics = eqns(abs0 = inpt),
                       thetas  = eqns(tvmtt = 1,
                                      tvntr = 1,
                                      tvf1 = 0.1),
                       omegas = eqns(etamtt = 0,
                                     etantr = 0,
                                     etaf1 = 0)){

  qpModel(params = eqns(MTT = exp(tvmtt + etamtt),
                        NTR = exp(tvntr + etantr),
                        LF1 = exp(tvf1 + etaf1)),
          theta = thetas,
          omega = omegas,
          output = eqns(MTT = MTT,
                        NTR = NTR,
                        LF1 = LF1,
                        etamtt = etamtt,
                        etantr = etantr,
                        etaf1 = etaf1),
          global = eqns(r = '// initialize dose amt and dose time arrays
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
X = 0.00001'),
          custom = eqns(r.main = '// TCAM variables; be sure to define NTR and MTT
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
 "END DO'),
          algebraic = algebraics)



}

#####  elimination models #####

# elimination
qp.elim = function(thetas = eqns(tvke = 1),
                   omegas = eqns(etake = 0)){
  qpModel(algebraic = eqns(elim = -ke*Acentral),
          params = eqns(ke = exp(tvke + etake)),
          theta = eqns(tvke = 1),
          omega = eqns(etake = 0),
          output = eqns(ke=ke,
                        etake=etake)
          )
}

# clearance
qp.elim.clearance = function(thetas = eqns(tvCL = 1,
                                           tvV = 1),
                             omegas = eqns(etaCL = 0,
                                           etaV = 0)){
  qpModel(algebraic = eqns(elim = -(CL/V)*Acentral),
             params = eqns(CL = exp(tvCL + etaCL),
                           V = exp(tvV + etaV)),
             theta = thetas,
             omega = omegas,
             output = eqns(CL=CL,
                           V=V,
                           etaCL=etaCL,
                           etaV=etaV)
             )
}

# Michaelis Menton clearance
qp.elim.MM = function(thetas = eqns(tvkmax = 1,
                                    tvkm50 = 1),
                      omegas = eqns(etakmax = 0,
                                    etakm50 = 0)){
  qpModel(algebraic = eqns(elim = -(kmax*Acentral)/(km50+Acentral)),
             params = eqns(kmax = exp(tvkmax + etakmax),
                           km50 = exp(tvkm50 + etakm50)),
             theta = thetas,
             omega = omegas,
             output = eqns(kmax=kmax,
                           km50=km50,
                           etakmax=etakmax,
                           etakm50=etakm50)
                     )
}

#####  compartmental models #####

# one cpt IV rate
qp.onecpt.iv.rate = function(thetas = eqns(tvke = 1,
                                      tvV = 1),
                        omegas = eqns(etake = 0,
                                      etaV = 0),
                        sigmas = eqns(add = 0,
                                      prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim),
          algebraic = eqns(abs0 = 0,
                           elim = ke*Acentral),
          params = eqns(ke = exp(tvke + etake),
                        V = exp(tvV + etaV)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        ke=ke,
                        V=V,
                        etake=etake,
                        etaV=etaV),
          custom = eqns(),
          global = eqns()
                       )
}

# one cpt IV cl
qp.onecpt.iv.cl = function(thetas = eqns(tvcl = 1,
                                           tvV = 1),
                             omegas = eqns(etacl = 0,
                                           etaV = 0),
                             sigmas = eqns(add = 0,
                                           prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim),
          algebraic = eqns(abs0 = 0,
                           elim = (cl/V)*Acentral),
          params = eqns(cl = exp(tvcl + etake),
                        V = exp(tvV + etaV)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        cl=cl,
                        V=V,
                        etacl=etacl,
                        etaV=etaV),
          custom = eqns(),
          global = eqns()
  )
}

# one cpt oral rate
qp.onecpt.oral.rate = function(thetas = eqns(tvke = 1,
                                        tvV = 1,
                                        tvka = 1),
                          omegas = eqns(etake = 0,
                                        etaV = 0,
                                        tvka = 0),
                          sigmas = eqns(add = 0,
                                        prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                      Acentral = abs0 - elim),
          algebraic = eqns(abs0 = ka*Aabs,
                           elim = ke*Acentral),
           params = eqns(ke = exp(tvke + etake),
                         ka = exp(tvka + etaka),
                         V = exp(tvV + etaV)),
           theta = thetas,
           omega = omegas,
           sigma = sigmas,
           observe = eqns(IPRED = Acentral/V,
                          DV = IPRED*(1+prop)+add),
           output = eqns(IPRED=IPRED,
                         DV=DV,
                         ke=ke,
                         V=V,
                         ka=ka,
                         etake=etake,
                         etaV=etaV)
    )
}


# one cpt oral cl
qp.onecpt.oral.cl = function(thetas = eqns(tvcl = 1,
                                             tvV = 1,
                                             tvka = 1),
                               omegas = eqns(etacl = 0,
                                             etaV = 0,
                                             etaka = 0),
                               sigmas = eqns(add = 0,
                                             prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                      Acentral = abs0 - elim),
          algebraic = eqns(abs0 = ka*Aabs,
                           elim = (cl/V)*Acentral),
          params = eqns(cl = exp(tvcl + etacl),
                        ka = exp(tvka + etaka),
                        V = exp(tvV + etaV)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        cl=cl,
                        V=V,
                        ka=ka,
                        etacl=etacl,
                        etaV=etaV,
                        etaka = etaka)
  )
}

# two cpt iv rate
qp.twocpt.iv.rate = function(thetas = eqns(tvke = 1,
                                           tvV = 1,
                                           tvk12 = 1,
                                           tvk21 = 1),
                             omegas = eqns(etake = 0,
                                           etaV = 0,
                                           etak12 = 0,
                                           etak21 = 0),
                             sigmas = eqns(add = 0,
                                           prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim - k12*Acentral + k21*Aperiph,
                      Aperiph = k12*Acentral - k21*Aperiph),
         algebraic = eqns(abs0 = 0,
                          elim = ke*Acentral),
         params = eqns(ke = exp(tvke + etake),
                       V = exp(tvV + etaV),
                       k12 = exp(tvk12 + etak12),
                       k21 = exp(tvk21 + etak21)),
         theta = thetas,
         omega = omegas,
         sigma = sigmas,
         observe = eqns(IPRED = Acentral/V,
                        DV = IPRED*(1+prop)+add),
         output = eqns(IPRED=IPRED,
                       DV=DV,
                       ke=ke,
                       V=V,
                       k12=k12,
                       k21=k21,
                       etake=etake,
                       etaV=etaV,
                       etak12=etak12,
                       etak21=etak21)
)
}

# two cpt iv cl
qp.twocpt.iv.cl = function(thetas = eqns(tvCL = 1,
                                         tvV = 1,
                                         tvQ = 1,
                                         tvV2 = 1),
                           omegas = eqns(etaCL = 0,
                                         etaV = 0,
                                         etaQ = 0,
                                         etaV2 = 0),
                           sigmas = eqns(add = 0,
                                         prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim - Q*(Acentral/V - Aperiph/V2),
                      Aperiph =Q*(Acentral/V - Aperiph/V2)),
          algebraic = eqns(abs0 = 0,
                           elim = (CL/V)*Acentral),
          params = eqns(CL = exp(tvCL + etaCL),
                        V = exp(tvV + etaV),
                        Q = exp(tvQ + etaQ),
                        V2 = exp(tvV2 + etaV2)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,CL=CL,
                        V=V,
                        Q=Q,
                        V2=V2,
                        etaCL=etaCL,
                        etaV=etaV,
                        etaQ=etaQ,
                        etaV2=etaV2)
  )
}

# two cpt oral rate
qp.twocpt.oral.rate = function(thetas = eqns(tvke = 1,
                                             tvV = 1,
                                             tvk12 = 1,
                                             tvk21 = 1,
                                             tvka = 1),
                               omegas = eqns(etake = 0,
                                             etaV = 0,
                                             etak12 = 0,
                                             etak21 = 0,
                                             etaka = 0),
                               sigmas = eqns(add = 0,
                                             prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                                        Acentral = abs0 - elim - k12*Acentral + k21*Aperiph,
                                        Aperiph = k12*Acentral - k21*Aperiph),
                            algebraic = eqns(abs0 = ka*Aabs,
                                             elim = ke*Acentral),
                            params = eqns(ke = exp(tvke + etake),
                                          V = exp(tvV + etaV),
                                          k12 = exp(tvk12 + etak12),
                                          k21 = exp(tvk21 + etak21),
                                          ka = exp(tvka + etaka)),
                            theta = thetas,
                            omega = omegas,
                            sigma = sigmas,
                            observe = eqns(IPRED = Acentral/V,
                                           DV = IPRED*(1+prop)+add),
                            output = eqns(IPRED=IPRED,
                                          DV=DV,
                                          ke=ke,
                                          V=V,k12=k12,
                                          k21=k21,
                                          ka=ka,
                                          etake=etake,
                                          etaV=etaV,
                                          etak12=etak12,
                                          etak21=etak21,
                                          etaka=etaka)
  )
}

# two cpt oral cl
qp.twocpt.oral.cl = function(thetas = eqns(tvCL = 1,
                                           tvV = 1,
                                           tvQ = 1,
                                           tvV2 = 1,
                                           tvka = 1),
                             omegas = eqns(etaCL = 0,
                                           etaV = 0,
                                           etaQ = 0,
                                           etaV2 = 0,
                                           etaka = 0),
                             sigmas = eqns(add = 0,
                                           prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                      Acentral = abs0 - elim - Q*(Acentral/V - Aperiph/V2),
                      Aperiph =Q*(Acentral/V - Aperiph/V2)),
          algebraic = eqns(abs0 = ka*Aabs,
                           elim = (CL/V)*Acentral),
          params = eqns(CL = exp(tvCL + etaCL),
                        V = exp(tvV + etaV),
                        Q = exp(tvQ + etaQ),
                        V2 = exp(tvV2 + etaV2),
                        ka = exp(tvka + etaka)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        CL=CL,V=V,
                        Q=Q,
                        V2=V2,
                        ka=ka,
                        etaCL=etaCL,
                        etaV=etaV,
                        etaQ=etaQ,
                        etaV2=etaV2,
                        etaka=etaka)
  )
}

# three cpt iv rate
qp.threecpt.iv.rate = function(thetas = eqns(tvke = 1,
                                           tvV = 1,
                                           tvk12 = 1,
                                           tvk21 = 1,
                                           tvk13 = 1,
                                           tvk31 = 1),
                             omegas = eqns(etake = 0,
                                           etaV = 0,
                                           etak12 = 0,
                                           etak21 = 0,
                                           etak13 = 0,
                                           etak31 = 0),
                             sigmas = eqns(add = 0,
                                           prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim - k12*Acentral + k21*Aperiph - k13*Acentral + k31*A2periph,
                      Aperiph = k12*Acentral - k21*Aperiph,
                      A2periph = k13*Acentral - k31*A2periph),
          algebraic = eqns(abs0 = 0,
                           elim = ke*Acentral),
          params = eqns(ke = exp(tvke + etake),
                        V = exp(tvV + etaV),
                        k12 = exp(tvk12 + etak12),
                        k21 = exp(tvk21 + etak21),
                        k13 = exp(tvk13 + etak13),
                        k31 = exp(tvk31 + etak31)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        ke=ke,
                        V=V,
                        k12=k12,
                        k21=k21,
                        k13=k13,
                        k31=k31,
                        etake=etake,
                        etaV=etaV,
                        etak12=etak12,
                        etak21=etak21)
  )
}

# three cpt iv cl
qp.threecpt.iv.cl = function(thetas = eqns(tvCL = 1,
                                         tvV = 1,
                                         tvQ = 1,
                                         tvV2 = 1,
                                         tvQ2 = 1,
                                         tvV3 = 1),
                           omegas = eqns(etaCL = 0,
                                         etaV = 0,
                                         etaQ = 0,
                                         etaV2 = 0,
                                         etaQ2 = 0,
                                         etaV3 = 0),
                           sigmas = eqns(add = 0,
                                         prop = 0)){
  qpModel(odes = eqns(Acentral = abs0 - elim - Q*(Acentral/V - Aperiph/V2) - Q2*(Acentral/V - A2periph/V3),
                      Aperiph = Q*(Acentral/V - Aperiph/V2),
                      A2periph = Q2*(Acentral/V - A2periph/V3)),
          algebraic = eqns(abs0 = 0,
                           elim = (CL/V)*Acentral),
          params = eqns(CL = exp(tvCL + etaCL),
                        V = exp(tvV + etaV),
                        Q = exp(tvQ + etaQ),
                        V2 = exp(tvV2 + etaV2),
                        Q2 = exp(tvQ2 + etaQ2),
                        V3 = exp(tvV3 + etaV3)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,CL=CL,
                        V = V,
                        Q = Q,
                        V2 = V2,
                        Q2 = Q2,
                        V3 = V3,
                        etaCL = etaCL,
                        etaV = etaV,
                        etaQ = etaQ,
                        etaV2 = etaV2)
  )
}

# three cpt oral rate
qp.threecpt.oral.rate = function(thetas = eqns(tvke = 1,
                                             tvV = 1,
                                             tvk12 = 1,
                                             tvk21 = 1,
                                             tvk13 = 1,
                                             tvk31 = 1,
                                             tvka = 1),
                               omegas = eqns(etake = 0,
                                             etaV = 0,
                                             etak12 = 0,
                                             etak21 = 0,
                                             etak13 = 0,
                                             etak31 = 0,
                                             etaka = 0),
                               sigmas = eqns(add = 0,
                                             prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                      Acentral = abs0 - elim - k12*Acentral + k21*Aperiph - k13*Acentral + k31*A2periph,
                      Aperiph = k12*Acentral - k21*Aperiph,
                      A2periph = k13*Acentral - k31*A2periph),
          algebraic = eqns(abs0 = ka*Aabs,
                           elim = ke*Acentral),
          params = eqns(ke = exp(tvke + etake),
                        V = exp(tvV + etaV),
                        k12 = exp(tvk12 + etak12),
                        k21 = exp(tvk21 + etak21),
                        k13 = exp(tvk13 + etak13),
                        k31 = exp(tvk31 + etak31),
                        ka = exp(tvka + etaka)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        ke=ke,
                        V=V,
                        k12=k12,
                        k21=k21,
                        k13 = k13,
                        k31 = k31,
                        ka=ka,
                        etake=etake,
                        etaV=etaV,
                        etak12=etak12,
                        etak21=etak21,
                        etak13 = etak13,
                        etak31 = etak31,
                        etaka=etaka)
  )
}

# three cpt oral clearance
qp.threecpt.oral.cl = function(thetas = eqns(tvCL = 1,
                                           tvV = 1,
                                           tvQ = 1,
                                           tvV2 = 1,
                                           tvQ2 = 1,
                                           tvV3 = 1,
                                           tvka = 1),
                             omegas = eqns(etaCL = 0,
                                           etaV = 0,
                                           etaQ = 0,
                                           etaV2 = 0,
                                           etaQ2 = 0,
                                           etaV3 = 0,
                                           etaka = 0),
                             sigmas = eqns(add = 0,
                                           prop = 0)){
  qpModel(odes = eqns(Aabs = -abs0,
                      Acentral = abs0 - elim - Q*(Acentral/V - Aperiph/V2) - Q2*(Acentral/V - A2periph/V3),
                      Aperiph = Q*(Acentral/V - Aperiph/V2),
                      A2periph = Q2*(Acentral/V - A2periph/V3)),
          algebraic = eqns(abs0 = ka*Aabs,
                           elim = (CL/V)*Acentral),
          params = eqns(CL = exp(tvCL + etaCL),
                        V = exp(tvV + etaV),
                        Q = exp(tvQ + etaQ),
                        V2 = exp(tvV2 + etaV2),
                        Q2 = exp(tvQ2 + etaQ2),
                        V3 = exp(tvV3 + etaV3),
                        ka = exp(tvka + etaka)),
          theta = thetas,
          omega = omegas,
          sigma = sigmas,
          observe = eqns(IPRED = Acentral/V,
                         DV = IPRED*(1+prop)+add),
          output = eqns(IPRED=IPRED,
                        DV=DV,
                        CL=CL,
                        V=V,
                        Q=Q,
                        V2=V2,
                        Q2=Q2,
                        V3=V3,
                        ka=ka,
                        etaCL=etaCL,
                        etaV=etaV,
                        etaQ=etaQ,
                        etaV2=etaV2,
                        etaka=etaka)
  )
}


##########   update functions ##########

# Update Theta values
theta = function(thetas = eqns()){
  qpModel(theta = thetas)
}

# Update Omega values
omega = function(omegas = eqns()){
  qpModel(omega = omegas)
}

# Update Sigma values
sigma = function(sigmas = eqns()){
  qpModel(sigma = sigmas)
}

# Update Parameter definitions
parameter = function(parameters = eqns()){
  qpModel(params = parameters)
}

odes = function(odes = eqns()){
  qpModel(odes = odes)
}
