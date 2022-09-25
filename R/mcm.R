#' Estimate and Test Inter-generational Social Mobility
#' Effect on an Outcome
#'
#' This function implements the mobility contrast model
#' designed for estimating and testing inter-generational
#' mobility effect on an outcome.
#'
#' @param formula an object of class "formula" (or one that can
#' be coerced to that class): a symbolic description of the model
#' to be fitted. A typical model used in studying social mobility
#' takes the form \code{response ~ origin*destination}, where
#' \code{respose} is the numeric response vector and \code{origin}
#' (\code{destination}) is a vector indicating the origin (destination).
#' The specification of \code{origin*destination} indicates the cross of
#' \code{origin} and \code{destination}, which is the same as \code{
#' origin + destination + origin:destination} where
#' \code{origin:destination} indicates the interaction of \code{origin}
#' and \code{destination}.
#'
#' @param data an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model. If not found in data,
#' the variables are taken from environment(formula),
#' typically the environment from which the function is called.
#' @param origin a character indicating the column name of origin.
#' @param destination a character indicating the column name of destination.
#' @param weights an optional vector of unit-level sampling weights to
#' be used in analysis.
#' Should be NULL or a numeric vector.
#' @param na.action a function which indicates what should
#' happen when the data contain NAs.The default is set by the
#' \code{na.action} setting in \code{options} and is
#' \code{na.fail} if that is unset.
#' @param family a character string, a function or the result of a call
#' to a family function describing the error distribution and link function
#' to be used in the model.
#' @param contrasts an optional list. The default is set as sum-to-zero
#' contrast.
#' @param gee logical. Should gee be used in estimating the model?
#' @param id a vector which identifies the clusters, which is required while
#' \code{gee} is used. The length of \code{id} should be the same as
#' the number of observations. Data are assumed to be sorted
#' so that observations on a cluster are contiguous rows for
#' all entities in the formula.
#' @param corstr a character string specifying
#' the correlation structure.
#' The following are permitted: \code{"independence"},
#' \code{"fixed"},
#' \code{"stat_M_dep"},
#' \code{"non_stat_M_dep"}, \code{"exchangeable"},
#' \code{"AR-M"} and \code{"unstructured"}.
#' @param displayresult logical. Should model results be displayed
#' after estimation. The default is \code{TRUE}.
#' @param \dots additional arguments to be passed to the function.
#'
#' @return A list containing:
#' \item{model}{Fitted generalized models of outcome on predictors.
#' See more on function \code{glm} in package \code{stats}.}
#' \item{origin_main}{Estimated main effects of origin.}
#' \item{destination_main}{Estimated main effects of destination.}
#' \item{mobility_estimates}{Estimated mobility effects.}
#' \item{mobility_se}{Standard errors of the estimated mobility effects.}
#' \item{mobility_sig}{Statistical significance of the the
#' estimated mobility effects.}
#'
#' @examples
#' library(MCM)
#' data('sim_moderate_het')
#' mcm(response ~ origin * destination, data = sim_moderate_het,
#'     origin = "origin",destination="destination")
#' @export

# mcm function used to estimate the mobility effect
mcm <- function(formula, data, weights=1, na.action=na.omit,
                origin,destination,family = gaussian(),
                contrasts = NULL,
                # contrasts = list(origin = "contr.sum",origin = "contr.sum"),
                gee = FALSE, # for gee extension
                id = NULL,
                corstr = "exchangeable",
                displayresult = TRUE,
                ...){
  op <- options(contrasts=c("contr.sum","contr.sum"), na.action = na.omit)
  on.exit(options(op))

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  fam <- family
  # parse formula
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # convert to factors
  mf[,origin] <- as.factor(mf[,origin])
  mf[,destination] <- as.factor(mf[,destination])

  mf[,"origin"] <- as.factor(mf[,origin])
  mf[,"destination"] <- as.factor(mf[,destination])

  mt <- attr(mf, "terms")

  if (!is.empty.model(mt)){
    x <- model.matrix(mt, mf, contrasts)
  }

  y <- model.response(mf, "numeric")

  # construct mobility variables
  # mobility status; 1=mobile, 0=nonmobile
  mf$mobility = as.factor(as.numeric(mf[,destination]!=mf[,origin]))
  # 1 = downward mobility; 2 = upward mobility
  mf = dplyr::mutate(mf,dir = as.factor(dplyr::case_when(as.numeric(mf$destination)==as.numeric(mf$origin) ~ 0,
                                                         as.numeric(mf$destination)<as.numeric(mf$origin) ~ 1,
                                                         as.numeric(mf$destination)>as.numeric(mf$origin) ~ 2)))
  # 1 = 1-step mobility; 2 = 2-step mobility
  mf$step <- 0
  mf = dplyr::mutate(mf,step = as.factor(dplyr::case_when(abs(as.numeric(mf$destination)-as.numeric(mf$origin))==1 ~ 1,
                                                          abs(as.numeric(mf$destination)-as.numeric(mf$origin))==2 ~ 2)))

  Orig = nlevels(mf$origin)
  Desti = nlevels(mf$destination)

  # formula <- update(formula, as.formula("~ . + origin*destination"))
  # estimate glm model
  if(gee==TRUE){
    temp6 = gee::gee(formula,
                     id = id,#id = get(id),
                     data = mf,
                     family = fam,
                     corstr = corstr)
  }else{
    # temp6 = glm(formula, mf, family = fam)
    temp6 = survey::svyglm(formula, design=survey::svydesign(ids=~1,
                                                             strata = NULL,
                                                             weights=~weights,
                                                 data=mf), family=fam)
  }
  # compute transformation matrix
  # trans.matrix = model.matrix(mt, mf, contrasts)
  # trans.matrix = trans.matrix[,stringr::str_detect(colnames(trans.matrix),":")]
  # trans.matrix = trans.matrix[!duplicated(trans.matrix),]

  twoway <- paste0("^",origin,"([0-9]*):",destination,"([0-9]*)$")

  mf <- mf[order(mf$destination),]
  mf <- mf[order(mf$origin),]

  trans.matrix = model.matrix(as.formula(paste0("~",origin,"*",destination)),mf)
  maino_trans.matrix = trans.matrix[,stringr::str_subset(colnames(trans.matrix), "^origin([0-9]*)$")]
  maino_trans.matrix = maino_trans.matrix[!duplicated(maino_trans.matrix),]
  maind_trans.matrix = trans.matrix[,stringr::str_subset(colnames(trans.matrix), "^destination([0-9]*)$")]
  maind_trans.matrix = maind_trans.matrix[!duplicated(maind_trans.matrix),]
  trans.matrix = trans.matrix[,stringr::str_subset(colnames(trans.matrix),twoway )]
  trans.matrix = trans.matrix[!duplicated(trans.matrix),]

  # all interaction estimates and se's
  if(gee==TRUE){
    ia_1 = temp6$coefficients[c(grep(paste0(":",destination), rownames(temp6$robust.variance)))]
    ia_2 = temp6$robust.variance[c(grep(paste0(":",destination), rownames(temp6$robust.variance))),c(grep(paste0(":",destination), rownames(temp6$robust.variance)))]
    maino_1 = temp6$coefficients[c(grep(paste0("^",origin,"[0-9]*$"), rownames(temp6$robust.variance)))]
    maino_2 = temp6$robust.variance[c(grep(paste0("^",origin,"[0-9]*$"), rownames(temp6$robust.variance))),c(grep(paste0("^",origin,"[0-9]*$"), rownames(temp6$robust.variance)))]
    maind_1 = temp6$coefficients[c(grep(paste0("^",destination,"[0-9]*$"), rownames(temp6$robust.variance)))]
    maind_2 = temp6$robust.variance[c(grep(paste0("^",destination,"[0-9]*$"), rownames(temp6$robust.variance))),c(grep(paste0("^",destination,"[0-9]*$"), rownames(temp6$robust.variance)))]
  }else{
    ia_1 = temp6$coefficients[c(grep(paste0(":",destination), rownames(vcov(temp6))))]
    ia_2 = vcov(temp6)[c(grep(paste0(":",destination), rownames(vcov(temp6)))),c(grep(paste0(":",destination), rownames(vcov(temp6))))]
    maino_1 = temp6$coefficients[c(grep(paste0("^",origin,"[0-9]*$"), rownames(vcov(temp6))))]
    maino_2 = vcov(temp6)[c(grep(paste0("^",origin,"[0-9]*$"), rownames(vcov(temp6)))),c(grep(paste0("^",origin,"[0-9]*$"), rownames(vcov(temp6))))]
    maind_1 = temp6$coefficients[c(grep(paste0("^",destination,"[0-9]*$"), rownames(vcov(temp6))))]
    maind_2 = vcov(temp6)[c(grep(paste0("^",destination,"[0-9]*$"), rownames(vcov(temp6)))),c(grep(paste0("^",destination,"[0-9]*$"), rownames(vcov(temp6))))]
  }

  # trans.matrix <- trans.matrix[,match(names(ia_1),colnames(trans.matrix))]
  iaesti = as.vector(trans.matrix%*%ia_1)
  iavcov = trans.matrix%*%ia_2%*%t(trans.matrix)
  maino_esti = as.vector(maino_trans.matrix%*%maino_1)
  maino_vcov = sqrt(diag(maino_trans.matrix%*%maino_2%*%t(maino_trans.matrix)))
  mtp_o = pt(-abs(maino_esti/maino_vcov), temp6$df.residual)*2 #p-values
  mtsig_o = rep('   ', Orig); mtsig_o[mtp_o<.05] = '*  '; mtsig_o[mtp_o<.01] = '** '; mtsig_o[mtp_o<.001] = '***'
  maino = data.frame(main_effects_origin=maino_esti,se=maino_vcov,sig=mtsig_o)
  rownames(maino) = paste0("Orig",1:nrow(maino))
  maind_esti = as.vector(maind_trans.matrix%*%maind_1)
  maind_vcov = sqrt(diag(maind_trans.matrix%*%maind_2%*%t(maind_trans.matrix)))
  mtp_d = pt(-abs(maind_esti/maind_vcov), temp6$df.residual)*2 #p-values
  mtsig_d = rep('   ', Desti); mtsig_d[mtp_d<.05] = '*  '; mtsig_d[mtp_o<.01] = '** '; mtsig_d[mtp_d<.001] = '***'
  maind = data.frame(main_effects_destination=maind_esti,se=maind_vcov,sig=mtsig_d)
  rownames(maind) = paste0("Desti",1:nrow(maind))

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
  sqrt(diag(m%*%(iavcov[(Desti*(i-1)+1):(Desti*i),(Desti*(i-1)+1):(Desti*i)])%*%t(m)))
})
mtse <- matrix(unlist(byrow_matrix_se),Orig,Desti,byrow = TRUE)

  # compute mobility effect significance
  mtp = pt(-abs(mtesti/mtse), temp6$df.residual)*2 #p-values
  mtp = matrix(mtp, Orig,Desti)

  mtsig = rep('   ', Orig*Desti); mtsig[mtp<.05] = '*  '; mtsig[mtp<.01] = '** '; mtsig[mtp<.001] = '***'
  mtsig = matrix(mtsig, Orig,Desti)

  mtesti_viz = matrix(sprintf("%.3f",round(mtesti,3)),Orig,Desti,byrow = FALSE)
  mtse_viz = matrix(sprintf("%.3f",round(mtse,3)),Orig, Desti, byrow = FALSE)
  mtsig_viz = mtsig
  diag(mtse)  ='-'; diag(mtesti)='-'; diag(mtp)   ='-'; diag(mtsig) ='-'
  diag(mtse_viz)  ='-----'; diag(mtesti_viz)='-----'; diag(mtsig_viz)   =''

  # display results
if(displayresult==TRUE){
  cat(
    # title
    "\n","      ",
    stringr::str_pad(paste0(" Desti ",1:Desti),10,"right"),
    # estimates, se, and significance
    sapply(1:Orig,function(i){
      c("\n",paste0("Orig ",i),
        stringr::str_pad(paste0(stringr::str_pad(mtesti_viz[i,],7,"left"," "),mtsig_viz[i,]),10,"right"),
        "\n","      ",
        stringr::str_pad(paste0("(",mtse_viz[i,],")"),10,"both")
      )
    }),sep = "    "
  )
}

  list(model = temp6,
       origin_main = maino,
       destination_main = maind,
       mobility_estimates = mtesti,
       mobility_se = mtse,
       mobility_sig = mtsig)
}
