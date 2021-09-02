
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
  ia_1 = broomExtra::tidy(model)$estimate[c(grep(twoway, broomExtra::tidy(model)$term ))]
  names(ia_1) <- broomExtra::tidy(model)$term[c(grep(twoway, broomExtra::tidy(model)$term ))]
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
  mtp = pt(-abs(mtesti/mtse), (nrow(model@frame)-nrow( broomExtra::tidy(model) ) + 1) )*2 #p-values
  mtp = matrix(mtp, Orig,Desti)

  mtsig = rep('   ', Orig*Desti); mtsig[mtp<.05] = '*  '; mtsig[mtp<.01] = '** '; mtsig[mtp<.001] = '***'
  mtsig = matrix(mtsig, Orig,Desti)

  diag(mtse)  ='-'; diag(mtesti)='-'; diag(mtp)   ='-'; diag(mtsig) ='-'

  # compute transformation matrix for three-way interactions -----
  # if age is continuous, there is no need to fill in the missed category
  threeway <- paste0("^",origin,"([0-9]*):",destination,"([0-9]*):",time,"$","|","^",time,":",origin,"([0-9]*):",destination,"([0-9]*)$")
  ia3_1 = broomExtra::tidy(model)$estimate[stringr::str_which(broomExtra::tidy(model)$term,threeway)]
  ia3_2 = vcov(model)[stringr::str_subset(broomExtra::tidy(model)$term,threeway),
                      stringr::str_subset(broomExtra::tidy(model)$term,threeway)]
  # broomExtra::tidy(model)[str_which(broomExtra::tidy(model)$term,threeway),]
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
  mtp = pt(-abs(mtesti3/mtse3), (nrow(model@frame)-nrow(broomExtra::tidy(model)) + 1) )*2 #p-values
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


