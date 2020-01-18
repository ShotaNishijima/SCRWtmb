#'
#' @import TMB
#'
NULL

#' Function for Estimating Correlated Random-Walk from Tracking Data
#'
#' @importFrom TMB MakeADFun
#' @importFrom TMB sdreport
#' @param data tracking data with the first column of longitude (easting) and latitude (northing)
#' @param mode whether different modes (migration and staying) is assumed (TRUE) or not (FALSE; default)
#' @param obs_error_type observation error distribution "normal" (default) or "t" (NOT work at present)
#' @param overdispersion whether overdispersion in binomial distribution is estimated (TRUE; may not walk) or not (FALSE; default)
#' @param minDF minimum df parameter in t-distribution
#' @param gamma_fix set fixed gamma value for each mode if neccessary (e.g., c(1,0))
#' @param theta_fix set fixed theta value for each mode if neccessary (e.g., c(0,0))
#' @param sigma_key able to assume the same observation error bewteen longitude and latitude by using c(0,0)
#' @param mixture whether using mixture of two modes
#' @param logit_p_init initial value for logit_p
#' @param phase the number of phase estimation when mode=TRUE
#' @encoding UTF-8
#' @export
scrw = function(data,
                mode = FALSE, # if TRUE, immigration mode and staying mode
                obs_error_type = "normal", # or "normal"
                overdispersion = FALSE,
                minDF = 2,
                gamma_fix = rep(-1,as.numeric(mode)+1), # please set c(NA, 0) if
                theta_fix = rep(-1,as.numeric(mode)+1),
                sigma_key = c(0,1), # 0:1 if assuming different observation error between longitude and latitude
                mixture = TRUE,
                logit_p_init = -0.5,
                phase = 1 # number of phase at mode=1
){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  if (class(data) != "matrix") data = as.matrix(data)
  y = as.matrix(stats::na.omit(data))
  timeSteps = 1:nrow(data)-1
  timeSteps = timeSteps[sapply(1:nrow(data), function(i) all(!is.na(data[i,])))]
  mode = as.numeric(mode)
  if (obs_error_type == "t") {
    obs_error_type_tmb = 0
  } else {
    obs_error_type_tmb = 1
  }
  if (length(unique(sigma_key))==1) {
    sigma_key = c(0,0)
  } else {
    sigma_key = c(0,1)
  }
  # make input data
  data_list = list(y=y,obs_t=timeSteps,timeSteps=max(timeSteps)+1,minDF=minDF,
                   mode = mode,gamma_fix = gamma_fix,theta_fix=theta_fix,
                   obs_error_type=obs_error_type_tmb,sigma_key=sigma_key,mixture = as.numeric(mixture))

  # set parameter initial values
  param = list()
  alpha = rep(0.8,mode+1)
  param$logit_alpha = log(alpha/(1-alpha))
  param$log_phi = 0
  if (mode==0) {
    gamma = 0.5
  } else {
    gamma = c(0.8,0.2)
  }
  param$logit_gamma = log(gamma)-log(1-gamma)
  param$trans_theta = rep(0,mode+1)
  param$log_sigma = rep(-1,sigma_key[2]+1)
  param$trans_rho = 0
  param$trans_df = ifelse(obs_error_type=="t",log(10-minDF),0)
  param$logit_p = rep(logit_p_init,max(timeSteps)+1)
  x_init = data
  x_init[is.na(x_init)] <- 0
  x_init <- x_init + 0.001
  param$x = x_init

  # set parameter(s) not to be estimated
  map = list()
  if (mode==0) {
    map$logit_p = rep(factor(NA),length(param$logit_p))
    map$logit_alpha = rep(factor(NA), length(param$logit_alpha))
    map$log_phi = factor(NA)
  }
  if (!isTRUE(overdispersion)) map$log_phi = factor(NA)
  map$logit_gamma = 0:(length(param$logit_gamma)-1)
  map$logit_gamma[gamma_fix != -1] <- NA
  map$logit_gamma <- factor(map$logit_gamma)
  map$trans_theta = 0:(length(param$trans_theta)-1)
  map$trans_theta[theta_fix != -1] <- NA
  map$trans_theta <- factor(map$trans_theta)
  map$trans_theta[gamma_fix==0] <- factor(NA)

  if (mode == 1 && phase >0) { # fix logit_p
    map$logit_p = rep(factor(NA),length(param$logit_p))
  }

  if (mode == 0) {
    f = TMB::MakeADFun(data=data_list, parameters=param, random = c("x"),
                  map = map, DLL="scrw")
  } else {
    f = TMB::MakeADFun(data=data_list, parameters=param, random = c("logit_p","x"),
                  map = map, DLL="scrw")
  }

  fit = stats::nlminb(f$par, f$fn, f$gr)
  p0_list = f$env$parList(fit$par)

  if (mode == 1 && phase >1) { # estimate only alpha and logit_p
    map2 = list()
    map2$log_phi = factor(NA)
    map2$logit_gamma = rep(factor(NA), length(param$logit_gamma))
    map2$trans_theta = rep(factor(NA), length(param$trans_theta))
    map2$log_sigma = rep(factor(NA), length(param$log_sigma))
    map2$trans_rho = factor(NA)
    map2$trans_df = factor(NA)
    # map2$logit_alpha = rep(factor(NA), length(param$logit_alpha))
    f = TMB::MakeADFun(data=data_list, parameters=p0_list, random = c("logit_p","x"),
                  map = map2, DLL="scrw")
    fit = stats::nlminb(f$par, f$fn, f$gr)
    p0_list = f$env$parList(fit$par)
  }

  if (mode == 1 && phase > 0) {
    map$logit_p = NULL
    f = TMB::MakeADFun(data=data_list, parameters=p0_list, random = c("logit_p","x"),
                  map = map, DLL="scrw")
    fit = stats::nlminb(f$par, f$fn, f$gr)
  }

  rep = TMB::sdreport(f, bias.correct = FALSE)
  p0_list = f$env$parList(fit$par)

  Res = list()
  Res$input <- arglist
  Res$opt <- fit
  Res$MakeADFun <- f
  Res$rep = rep
  Res$par_list = p0_list

  x = rep$par.random[names(rep$par.random) == "x"]
  x = matrix(x, ncol = ncol(param$x),nrow = nrow(param$x))
  Res$x <- x

  if (mode){
    p = rep$par.random[names(rep$par.random) == "logit_p"]
    p = 1/(1+exp(-p))
  } else {
    p = NULL
  }
  Res$p <- p

  if (mode) {
    Res$alpha = alpha <- rep$value[names(rep$value) == "alpha"]
    Res$phi = rep$value[names(rep$value) == "phi"]
    Res$mean_p = as.numeric(-(1-alpha[2])/(alpha[1]+alpha[2]-2))
  } else {
    Res$alpha <- Res$phi <- Res$mean_p <- NULL
  }

  Res$gamma = rep$value[names(rep$value) == "gamma"]
  Res$theta = rep$value[names(rep$value) == "theta"]
  Res$sigma = rep$value[names(rep$value) == "sigma"]
  Res$rho = rep$value[names(rep$value) == "rho"]
  if (obs_error_type == "t") {
    Res$df <- rep$value[names(rep$value) == "df"]
    Res$omega <- NULL
  } else {
    Res$df <- NULL
    Res$omega <- rep$value[names(rep$value) == "omega"]
  }

  NN = nrow(data)*2
  Res$loglik <- loglik <- -fit$objective
  Res$k <- k <- length(fit$par)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(NN-k-1)
  Res$BIC <- -2*loglik+k*log(NN)

  return(Res)
}


#' Function for activating cpp file for TMB
#'
#' @importFrom TMB compile
#' @importFrom TMB dynlib
#' @param TmbFile Cpp file name
#' @param CppDir directory having the cpp file
#' @param RunDir directory for running the code (default: working directory)
#' @param overwrite whether the cpp file is overwritten (TRUE) or not (FALSE; default)
#'
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' use_scrw_tmb()
#' }
#'
#' @export

use_scrw_tmb <- function(TmbFile = "scrw",
                         CppDir = system.file("executable",package="SCRWtmb"),
                         RunDir = getwd(),
                         overwrite = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Please install TMB package!")
  }
  file.copy( from=paste0(CppDir,"/",TmbFile,".cpp"), to=paste0(RunDir,"/",TmbFile,".cpp"), overwrite=overwrite)
  TMB::compile( paste0(TmbFile,".cpp") )
  dyn.load(TMB::dynlib(TmbFile))
}
