#' Estimate and Test Inter-generational Mobility Effect
#' with Longitudinal Data
#'
#' This function fits a multilevel mobility contrast model
#' to estimate and test inter-generational mobility effect
#' on an outcome in longitudinal data.
#'
#'
#' @param formula Inherit the function form from \code{lme4}
#' package. It is a two-sided linear formula
#' object describing both the fixed-effects and
#' random-effects part of the model, with the response
#' on the left of a \code{~} operator and the terms,
#' separated by + operators, on the right.
#' Random-effects terms are distinguished by vertical
#' bars (\code{|}) separating expressions for design
#' matrices from grouping factors.
#' Two vertical bars (\code{||}) can be used to specify
#' multiple uncorrelated random effects for the same
#' grouping variable. (Because of the way it is implemented,
#' the \code{||}-syntax works only for design matrices
#' containing numeric (continuous) predictors;
#' to fit models with independent categorical effects,
#' see dummy or the lmer_alt function from the \code{afex}
#' package.) A typical model used in studying social mobility with
#' longitudinal data
#' takes the form \code{response ~ origin*destination + | id}, where
#' \code{respose} is the numeric response vector and \code{origin}
#' (\code{destination}) is a vector indicating the origin (destination).
#' The specification of \code{origin*destination} indicates the cross of
#' \code{origin} and \code{destination}, which is the same as \code{
#' origin + destination + origin:destination} where
#' \code{origin:destination} indicates the interaction of \code{origin}
#' and \code{destination}. \code{id} is a identifier for the clusters.
#' @param data an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model. If not found in data,
#' the variables are taken from environment(formula),
#' typically the environment from which the function is called.
#' @param REML logical. Should the estimates be chosen be optimize the
#' restricted log-likelihood (REML) criterial (as opposed to the
#' log-likelihood)?
#' @param origin a character indicating the column name of origin.
#' @param destination a character indicating the column name of destination.
#' @param time a character indicating the time when individual was observed
#' @param control Inherit from \code{lme4} package. It is a list (of correct
#' class, resulting from lmerControl() or glmerControl() respectively)
#' containing control parameters, including the nonlinear optimizer to
#' be used and parameters to be passed through to the nonlinear optimizer,
#' see the \code{lmerControl} documentation in \code{lme4} package for details.
#' @param start Inherit from \code{lme4} package. It is a named list of
#' starting values for the parameters in the model.
#' @param verbose Inherit from \code{lme4} package. It is an integer scalar.
#' If > 0 verbose output is generated during the optimization of the parameter
#' estimates. If > 1 verbose output is generated during the individual
#' penalized iteratively reweighted least squares (PIRLS) steps.
#' @param subset optional expression selecting the subset of the rows of data
#' to fit the model.
#' @param weights an optional vector of ‘prior weights’ to
#' be used in the fitting process.
#' Should be NULL or a numeric vector.
#' @param na.action a function which indicates what should
#' happen when the data contain NAs.The default is set by the
#' \code{na.action} setting in \code{options} and is
#' \code{na.fail} if that is unset.
#' @param offset Inherit from \code{lme4} package. This can be used
#' to specify an a priori known component to be included in the linear
#' predictor during fitting. This should be NULL or a numeric vector
#' of length equal to the number of cases. One or more offset
#' terms can be included in the formula instead or as well,
#' and if more than one is specified their sum is used.
#' @param contrasts an optional list. The default is set as sum-to-zero
#' contrast.
#' @param devFunOnly logical - return only the deviance evaluation function.
#' @param displayresult logical. Should model results be displayed
#' after estimation. The default is \code{TRUE}.
#' @param \dots additional arguments to be passed to the function.
#'
#' @return A list containing:
#' \item{model}{Fitted generalized models of outcome on predictors.
#' See more on function \code{glm} in package \code{stats}.}
#' \item{estimates}{Estimated mobility effects.}
#' \item{se}{Standard errors of the estimated mobility effects.}
#' \item{significance}{Statistical significance of the the
#' estimated mobility effects.}
#' \item{esti_3way}{Estimated mobility effects conditional
#' on specific age.}
#' \item{se_3way}{Standard errors of the estimated mobility
#' effects conditional specific age.}
#' \item{sig_3way}{Statistical significance of the the
#' estimated mobility effects conditional on age.}
#'
#'
#' @examples
#' library(MCM)
#' library(lme4)
#' data("sim_datlmer")
#' fit_mcm_lmer <- mcm_lmer(yij ~ origin*destination*age +
#'                            (1|id), data = sim_datlmer,
#'                          origin = "origin",
#'                          destination = "destination",
#'                          time = "age")
#'
#'




mcm_lmer <- function(formula, data = NULL,
                     REML = TRUE,
                     control = lme4::lmerControl(),
                     start = NULL, verbose = 0L,
                     subset, weights, na.action, offset,
                     contrasts = NULL, devFunOnly = FALSE,
                     origin = NULL,destination = NULL,time = NULL,
                     displayresult = TRUE,
                     ...){
  # data <- as.data.frame(data)
  # data$origin <- data[,"origin"]
  # data$destination <- data[,"destination"]

  # contrast in estimating the model
  op <- options(contrasts=c("contr.sum","contr.sum"), na.action = na.omit)
  on.exit(options(op))

  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate. = TRUE)
    control <- do.call(lme4::lmerControl, control)
  }
  mc$control <- control
  mc[[1]] <- quote(lme4::lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  # added the interaction
  # mcout$formula <- update(mcout$formula,as.formula("~origin*destination +."))

  lmod$formula <- NULL

  devfun <- do.call(mkLmerDevfun, c(lmod, list(start = start,
                                               verbose = verbose, control = control)))
  if (devFunOnly)
    return(devfun)
  if (identical(control$optimizer, "none"))
    stop("deprecated use of optimizer=='none'; use NULL instead")
  opt <- if (length(control$optimizer) == 0) {
    s <- getStart(start, environment(devfun)$pp)
    list(par = s, fval = devfun(s), conv = 1000, message = "no optimization")
  }
  else {
    optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
                 boundary.tol = control$boundary.tol, control = control$optCtrl,
                 verbose = verbose, start = start, calc.derivs = control$calc.derivs,
                 use.last.params = control$use.last.params)
  }
  cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
                  lbound = environment(devfun)$lower)
  model <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
           mc = mcout, lme4conv = cc)


  # get the mobility effect
  data <- model@frame
  data[,origin] <- as.factor(data[,origin])
  data[,destination] <- as.factor(data[,destination])

  data[,"origin"] <- as.factor(data[,origin])
  data[,"destination"] <- as.factor(data[,destination])

  Orig = nlevels(as.factor(data$origin))
  Desti = nlevels(as.factor(data$destination))

  # compute transformation matrix for two-way interactions -----
  twoway <- paste0("^",origin,"([0-9]*):",destination,"([0-9]*)$")

  data <- data[order(data$destination),]
  data <- data[order(data$origin),]

  # data <- dplyr::arrange(data,origin,destination)
  trans.matrix = model.matrix(as.formula(paste0("~",origin,"*",destination)),data)
  # trans.matrix = model.matrix(as.formula(paste0("~",origin,"*",destination)),data[order(data[,c(origin)],data[,c(destination)]),])
  trans.matrix = trans.matrix[,stringr::str_subset(colnames(trans.matrix),twoway )]
  # trans.matrix = trans.matrix[,stringr::str_detect(colnames(trans.matrix),":")]
  trans.matrix = trans.matrix[!duplicated(trans.matrix),]

  # all interaction estimates and se's
  coefficienttable <- data.frame(term = parameters::model_parameters(model)$Parameter,
                                 estimate = parameters::model_parameters(model)$Coefficient)
  ia_1 = coefficienttable$estimate[c(grep(twoway, coefficienttable$term ))]
  names(ia_1) <- coefficienttable$term[c(grep(twoway, coefficienttable$term ))]
  ia_2 = vcov(model)[c(grep(twoway,
                            rownames(vcov(model)))),c(grep(twoway, rownames(vcov(model))))]

  # trans.matrix <- trans.matrix[,match(names(ia_1),colnames(trans.matrix))]
  trans.matrix <- trans.matrix[,match(names(ia_1),colnames(trans.matrix))]
  iaesti = as.vector(trans.matrix%*%ia_1)
  iavcov = trans.matrix%*%ia_2%*%t(trans.matrix)

  # mobility contrast and get the mobility effect estimates and SEs
  byrow_matrix_esti <- lapply(1:Orig,function(i){
    m <- matrix(0,Orig,Desti)
    m[,i] <- 1;m[i,] <- NA
    m <- diag(Desti)-m
    t(m%*%iaesti[(Desti*(i-1)+1):(Desti*i)])
  })
  mtesti <- matrix(unlist(byrow_matrix_esti),Orig,Desti,byrow = TRUE)

  byrow_matrix_se <- lapply(1:Orig,function(i){
    m <- matrix(0,Orig,Desti)
    m[,i] <- 1;m[i,] <- NA
    m <- diag(Desti)-m
    tempiavcov <- iavcov[(Desti*(i-1)+1):(Desti*i),(Desti*(i-1)+1):(Desti*i)]
    sqrt(diag(as.matrix(m%*%(tempiavcov)%*%t(m))))
  })
  mtse <- matrix(unlist(byrow_matrix_se),Orig,Desti,byrow = TRUE)

  # compute mobility effect significance
  # ask this residual things
  mtp = pt(-abs(mtesti/mtse), (nrow(model@frame)-nrow( coefficienttable ) + 1) )*2 #p-values
  mtp = matrix(mtp, Orig,Desti)

  mtsig = rep('   ', Orig*Desti); mtsig[mtp<.05] = '*  '; mtsig[mtp<.01] = '** '; mtsig[mtp<.001] = '***'
  mtsig = matrix(mtsig, Orig,Desti)

  diag(mtse)  ='-'; diag(mtesti)='-'; diag(mtp)   ='-'; diag(mtsig) ='-'

  # compute transformation matrix for three-way interactions -----
  # if age is continuous, there is no need to fill in the missed category
  threeway <- paste0("^",origin,"([0-9]*):",destination,"([0-9]*):",time,"$","|","^",time,":",origin,"([0-9]*):",destination,"([0-9]*)$")
  ia3_1 = coefficienttable$estimate[stringr::str_which(coefficienttable$term,threeway)]
  ia3_2 = vcov(model)[stringr::str_subset(coefficienttable$term,threeway),
                      stringr::str_subset(coefficienttable$term,threeway)]
  # parameters::Coefficient(model)[str_which(parameters::Coefficient(model)$term,threeway),]
  iaesti3 = as.vector(trans.matrix%*%ia3_1)
  iavcov3 = trans.matrix%*%ia3_2%*%t(trans.matrix)

  # mobility contrast and get the mobility effect estimates and SEs
  byrow_matrix_esti <- lapply(1:Orig,function(i){
    m <- matrix(0,Orig,Desti)
    m[,i] <- 1;m[i,] <- NA
    m <- diag(Desti)-m
    t(m%*%iaesti3[(Desti*(i-1)+1):(Desti*i)])
  })
  mtesti3 <- matrix(unlist(byrow_matrix_esti),Orig,Desti,byrow = TRUE)

  byrow_matrix_se <- lapply(1:Orig,function(i){
    m <- matrix(0,Orig,Desti)
    m[,i] <- 1;m[i,] <- NA
    m <- diag(Desti)-m
    tempiavcov <- iavcov[(Desti*(i-1)+1):(Desti*i),(Desti*(i-1)+1):(Desti*i)]
    sqrt(diag(as.matrix(m%*%(tempiavcov)%*%t(m))))
  })
  mtse3 <- matrix(unlist(byrow_matrix_se),Orig,Desti,byrow = TRUE)

  # compute mobility effect significance
  # ask this residual things
  mtp = pt(-abs(mtesti3/mtse3), (nrow(model@frame)-nrow(coefficienttable) + 1) )*2 #p-values
  mtp3 = mtp = matrix(mtp, Orig,Desti)

  mtsig = rep('   ', Orig*Desti); mtsig[mtp<.05] = '*  '; mtsig[mtp<.01] = '** '; mtsig[mtp<.001] = '***'
  mtsig3 = matrix(mtsig, Orig,Desti)

  mtesti3_viz = matrix(sprintf("%.3f",round(mtesti3,3)),Orig,Desti,byrow = FALSE)
  mtse3_viz = matrix(sprintf("%.3f",round(mtse3,3)),Orig, Desti, byrow = FALSE)
  mtsig3_viz = mtsig3

  # mtesti3_viz = round(mtesti3,3);mtse3_viz = round(mtse3,3);mtsig3_viz = mtsig3
  # diag(mtse3_viz)  ='—'; diag(mtesti3_viz)='—'; diag(mtsig3_viz) ='—'
  # diag(mtse3_viz)  ='—————'; diag(mtesti3_viz)='—————'; diag(mtsig3_viz)   =''
  diag(mtse3_viz)  ='-----'; diag(mtesti3_viz)='-----'; diag(mtsig3_viz)   =''

  # display results
  if(displayresult==TRUE){
  cat(
    # title
    "","Conditional on age, the mobility effect under effect coding is:",
    # column name:
    "\n","      ",
    stringr::str_pad(paste0(" Desti ",1:Desti),10,"right"),
    # estimates, se, and significance
    sapply(1:Orig,function(i){
      c("\n",paste0("Orig ",i),
        stringr::str_pad(paste0(stringr::str_pad(mtesti3_viz[i,],7,"left"," "),mtsig3_viz[i,]),10,"right"),
        "\n","      ",
        stringr::str_pad(paste0("(",mtse3_viz[i,],")"),10,"both")
      )
    }),sep = "    "
  )
}
  # cat(
  #   # title
  #   "","conditional on age, the mobility effect is:",
  #   # column name:
  #   "\n","      ",
  #   stringr::str_pad(paste0("Desti ",1:Desti),10,"right"),
  #   # estimates, se, and significance
  #   sapply(1:Orig,function(i){
  #     c("\n",paste0("Orig ",i),
  #       stringr::str_pad(paste0(stringr::str_sub(mtesti3_viz[i,],1,6),mtsig3_viz[i,]),10,"right"),
  #       "\n","      ",
  #       stringr::str_pad(paste0("(",stringr::str_sub(mtse3_viz[i,],1,6),")"),10,"right")
  #     )
  #   }),sep = "    "
  # )

  diag(mtse3)  ='-'; diag(mtesti3)='-'; diag(mtp3)   ='-'; diag(mtsig3) ='-'
  list(model = model,
       estimates = mtesti,
       se = mtse,
       significance = mtsig,
       esti_3way = mtesti3,
       se_3way = mtse3,
       sig_3way = mtsig3)
}


