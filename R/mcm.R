
# mcm function used to estimate the mobility effect
mcm <- function(formula, data, weights, na.action=na.omit,
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

  fam <- family
  if (is.character(family))
    fam <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    fam <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
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
                     data = fm,
                     family = fam,
                     corstr = corstr)
  }else{
    temp6 = glm(formula, mf, family = fam)
  }
  # compute transformation matrix
  # trans.matrix = model.matrix(mt, mf, contrasts)
  # trans.matrix = trans.matrix[,stringr::str_detect(colnames(trans.matrix),":")]
  # trans.matrix = trans.matrix[!duplicated(trans.matrix),]

  twoway <- paste0("^",origin,"([0-9]*):",destination,"([0-9]*)$")

  mf <- mf[order(mf$destination),]
  mf <- mf[order(mf$origin),]

  trans.matrix = model.matrix(as.formula(paste0("~",origin,"*",destination)),mf)
  trans.matrix = trans.matrix[,stringr::str_subset(colnames(trans.matrix),twoway )]
  trans.matrix = trans.matrix[!duplicated(trans.matrix),]

  # all interaction estimates and se's
  if(gee==TRUE){
    ia_1 = temp6$coefficients[c(grep(paste0(":",destination), rownames(temp6$robust.variance)))]
    ia_2 = temp6$robust.variance[c(grep(paste0(":",destination), rownames(temp6$robust.variance))),c(grep(paste0(":",destination), rownames(temp6$robust.variance)))]
  }else{
    ia_1 = temp6$coefficients[c(grep(paste0(":",destination), rownames(vcov(temp6))))]
    ia_2 = vcov(temp6)[c(grep(paste0(":",destination), rownames(vcov(temp6)))),c(grep(paste0(":",destination), rownames(vcov(temp6))))]
  }

  # trans.matrix <- trans.matrix[,match(names(ia_1),colnames(trans.matrix))]
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
       estimates = mtesti,
       se = mtse,
       significance = mtsig)
}
