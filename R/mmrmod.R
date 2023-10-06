#' MMRM analysis without assuming orthogonality property
#'
#' \code{mmrmod} provides inference results for the mixed models for repeated
#' measures (MMRM) analysis without assuming orthogonality property between
#' fixed effect and variance-covariance parameters.
#'
#' @param data a data frame that may include \code{outcome}, \code{group},
#'   \code{time}, \code{id}, and specified covariate variables.
#' @param outcome a name of outcome (dependent) variable included in
#'   \code{data}.
#' @param group  a name of treatment group variable included in \code{data}.
#' @param time a name of time variable for repeated measurements included
#'   in \code{data}. The default is \code{NULL}.
#' @param subject a name of subject id variable for repeated measurements
#'   included in \code{data}. The default is \code{NULL}.
#' @param covariate a character vector with names of covariate variables
#'   included in \code{data}. The default is \code{NULL}.
#' @param covfactor an integer vector including nominal variable indicators for
#'   covariate variables. Nominal variable: \code{1}, continuous variable:
#'   \code{0}. The default is \code{NULL}.
#' @param type specify the covariance structure type from \code{c("UN", "CS",
#'   "AR(1)")}. The default is \code{"UN"}.
#' @param hetero a logical value specifying whether or not to use a unequal
#'   variance model among groups.
#'   Unequal variance: \code{TRUE}, equal variance: \code{FALSE}.
#'   Default is \code{TRUE}.
#' @param robust a logical value specifying whether or not to use robust
#'   inference with sandwich-type standard error.
#'   Robust inference: \code{TRUE}, naive inference: \code{FALSE}.
#'   Default is \code{TRUE}.
#' @param ssadjust an integer specifying whether small sample adjustment is
#'   performed. No adjustment: \code{0}, adjustment method 1: \code{1},
#'   adjustment method 2: \code{2}. For details on each method, see
#'   Maruo et al. (submitting). Default is \code{1}.
#' @param reml a logical value specifying whether to use the REML inference.
#'   If orthogonality is not assumed, the REML inference cannot be used, so
#'   if \code{REML=TRUE} is specified, only an \code{mmrm} object by
#'   the \code{mmrm} package is returned. Default is \code{FALSE}.
#' @param intcov a logical value specifying whether to include interactions
#'   between \code{time} and \code{covariate}. Default is \code{TRUE}.
#' @param conf.level a numeric value of the confidence level for the
#'   confidence intervals. Default is 0.95.
#' @param ... additional arguments to be passed to the \code{mmrm} function.
#'   For example, it can identify an optimizer in parameter estimation.
#'   Note that the Kenward-Roger method cannot be used in degree-of-freedom
#'   estimation unless \code{REML=TRUE} is specified.
#'   See \code{\link{mmrm}} function for details.
#'
#' @details The inference without assuming orthogonality property means that the
#' off-diagonal ("od") blocks of the variance-covariance matrix of the
#' parameters, i.e., the covariance between the fixed effects parameter and
#' the covariance parameter, are not zero.
#'
#' The usual linear mixed model program packages make inferences on the
#' assumption that the off-diagonal block is a zero matrix from the properties
#' of the normal distribution, however, this assumption is incorrect when
#' the normality of the error does not hold and/or the missing values are not
#' from MCAR mechanism.
#'
#' The function uses the \code{mmrm} function from the \code{mmrm} package
#' for parameter estimation.
#'
#' @return a list that includes following components
#' for the MMRM analysis:
#' \describe{
#'   \item{\code{marg.mean}}{a data frame including inference results of
#'         marginal mean for each level of \code{time} and \code{group}.}
#'   \item{\code{diff.mean}}{a data frame including inference results of
#'         difference of marginal mean between \code{group} at each \code{time}.
#'         }
#'   \item{\code{Vtheta}}{a variance-covariance matrix for MLE of
#'         the vector of all parameters.}
#'   \item{\code{mmrmObject}}{an "\code{mmrm}" object by \code{mmrm} package. }
#' }

#'
#' @references \itemize{
#'   \item Maruo K, Ishii R, Yamaguchi Y, Doi M, Gosho M (2020). A note on the
#'         bias of standard errors when orthogonality of mean and variance
#'         parameters is not satisfied in the mixed model for repeated measures
#'         analysis. Statistics in Medicine, 39, 1264-1274,
#'         \doi{10.1002/sim.8474}.\cr
#'   \item Maruo K, Ishii R, Yamaguchi Y, Gosho M (submitting). Small sample
#'         adjustment for inference without assuming orthogonality in MMRM
#'         analysis.
#'   \item Sabanes Bove D, Dedic J, Kelkhoff D, Kunzmann K, Lang B, Li L, Wang Y
#'         (2022). mmrm: Mixed models for repeated measures. R package version
#'         0.2.2, \url{https://CRAN.R-project.org/package=mmrm}.
#' }
#'
#' @seealso \code{\link{mmrm}}
#'
#' @examples
#' library(mmrm)
#' # fev_data from mmrm package
#' data(fev_data)
#'
#' res <- mmrmod(data = fev_data, outcome = FEV1, group = ARMCD, time = AVISIT,
#'               subject = USUBJID, covariate = c("RACE", "SEX", "FEV1_BL"),
#'               covfactor = c(1, 1, 0))
#' res$marg.mean
#' res$diff.mean
#'
#' @importFrom stats as.formula model.matrix pnorm pt qnorm qt ftable coef xtabs
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom mmrm mmrm mmrm_control VarCorr
#'
#'
#' @export

mmrmod <- function(data, outcome, group, time, subject, covariate = NULL,
                   covfactor = NULL, type = "UN", hetero = TRUE, robust = TRUE,
                   ssadjust = 1, reml = FALSE, intcov = TRUE, conf.level = 0.95,
                   ...){
  od <- TRUE

  data00 <- as.data.frame(data)
  data <- data.frame(y = data00[, deparse(substitute(outcome))],
                     groupc = data00[, deparse(substitute(group))],
                     timec = data00[, deparse(substitute(time))],
                     id = as.factor(data00[, deparse(substitute(subject))]))
  nc <- 0
  if (!is.null(covariate)) {
    nc <- length(covariate)
    for (i in 1:nc){
      if (covfactor[i] == 1){
        data[, covariate[i]] <- as.factor(data00[, covariate[i]])
      } else {
        data[, covariate[i]] <- data00[, covariate[i]]

      }
    }
  }
  tl <- names(table(data$timec))
  nt <- length(tl)
  gl <- names(table(data$groupc))
  ng <- length(gl)

  data$group <- NA
  for (i in 1:ng){
    data$group[data$groupc == gl[i]] <- i
  }
  data$time <- NA
  for (i in 1:nt){
    data$time[data$timec == tl[i]] <- i
  }
  data$time <- as.factor(data$time)
  data$group <- as.factor(data$group)
  data0 <<- data
  form0 <- "y ~ group * time"
  if (!is.null(covariate)) {
    if (intcov) {
      form0 <- paste(form0, "+ time * (")
    } else {
      form0 <- paste(form0, "+")
    }
    for (i in 1:nc) {
      if (i == 1) {
        form0 <- paste(form0, covariate[i])
      } else {
        form0 <- paste(form0, "+", covariate[i])
      }
    }
    if (intcov) {
      form0 <- paste(form0, ")")
    }

    formc <- paste("y ~ ", paste(covariate, collapse = " + "))
    formc <- as.formula(formc)
  }
  formula0 <- as.formula(form0)
  if (hetero){
    if (type == "UN"){
      form0 <- paste(form0,
                     "+ us(time | group / id)")
    }
    if (type == "CS"){
      form0 <- paste(form0,
                     "+ cs(time | group / id)")
    }
    if (type == "AR(1)"){
      form0 <- paste(form0,
                     "+ ar1(time | group / id)")
    }
  } else {
    if (type == "UN"){
      form0 <- paste(form0, "+ us(time | id)")
    }
    if (type == "CS"){
      form0 <- paste(form0, "+ cs(time | id)")
    }
    if (type == "AR(1)"){
      form0 <- paste(form0, "+ ar1(time | id)")
    }
  }
  formula <- as.formula(form0)
  fit <- mmrm(formula = formula, data = data, reml = reml, ...)
  if (reml) {
    res <- fit
    warning("If reml = TRUE is specified, standard error without assuming
            orthogonality propaty cannot be estimated. Only mmrm object is
            returned.")
  } else {
    if (fit$method %in% c("Kenward-Roger", "Kenward-Roger-Linear")){
      stop("Kenward-Roger method cannot be used unless REML = TRUE is
           specified.")
    }
    if (hetero){
      V <- VarCorr(fit)
    } else {
      V <- list()
      V[[1]] <- VarCorr(fit)
    }
    beta <- coef(fit)
    nb <- length(beta)

    if (hetero){
      ngs <- ng
    } else {
      ngs <- 1
    }

    alp <- list()
    alp0 <- c()
    for (i in 1:ngs){
      if (type == "UN") {
        alp[[i]] <- V[[i]][upper.tri(V[[i]], diag = TRUE)]
      }
      if (type == "CS") {
        alp[[i]] <- c(V[[i]][1, 1] - V[[i]][1, 2], V[[i]][1, 2])
      }
      if (type == "AR(1)") {
        alp[[i]] <- c(V[[i]][1, 1], V[[i]][1, 2] / V[[i]][1, 1])
      }
      alp0 <- c(alp0, alp[[i]])
    }
    ns <- length(alp0)

    if (hetero){
      nsg <- ns / ng
    } else {
      nsg <- ns
    }

    tr <- function(X) sum(diag(X))
    Dec2Bin <- function(num, digit=0){
      if(num <= 0 && digit <= 0){
        return(NULL)
      }else{
        return(append(Recall(num %/% 2, digit - 1), num %% 2))
      }
    }

    evars <- as.character(formula)[3]
    casenames <- names(table(data$id))
    msflg <- table(data$id, is.na(data$y))[, 1]
    N <- sum(msflg != 0)
    omis <- names(msflg)[msflg == 0]
    if (length(omis) > 0L){
      for (i in 1:length(omis)){
        data <- data[data$id != omis[i], ]
      }
    }

    data0 <- data[!duplicated(data$id),]
    data0$y <- NULL
    data01 <- data0
    data1 <- c()
    for (i in 1:nt){
      data01$time <- i
      data1 <- rbind(data1,data01)
    }
    data1$time <- as.factor(data1$time)
    data2 <- data[, c("id", "time", "y")]
    data <- merge(data2, data1, by = c("id", "time"), all = T)
    data01 <- data[!duplicated(data$id),]
    grp0 <- data01$group

    Ncc <- sum(msflg == nt)
    flg.na <- !is.na(data$y)
    n.data <- sum(flg.na)
    options(na.action = "na.pass")
    X <- model.matrix(formula0, data = data)
    options(na.action = "na.omit")
    if (!is.null(covariate)){
      options(na.action = "na.pass")
      Xcov <- as.matrix(model.matrix(formc, data = data01))[, -1]
      options(na.action = "na.omit")
      covmean <- colMeans(Xcov)
    } else {
      covmean <- NULL
    }
    y <- data$y[flg.na]
    n.na <- length(flg.na)-n.data
    adj.prm <- c(N, Ncc, n.data, n.na)
    Xb <- list()
    idt <- cbind((1:N) %x% (numeric(nt) + 1), (numeric(N) + 1) %x% (1:nt))
    for (j in 1:nb){
      dum <- as.data.frame(cbind(idt, X[, j]))
      names(dum) <- c("a", "b", "c")
      Xb[[j]] <- as.matrix(ftable(xtabs(c ~ ., dum)))
    }
    dimnames(X) <- NULL
    nti <- table(data$id[!is.na(data$y)])
    ntic <- cumsum(nti)
    ntic <- c(0, ntic)
    X <- X[flg.na, ]
    re <- y - X %*% beta
    dum <- as.data.frame(cbind(idt[flg.na, ], re))
    names(dum) <- c("a", "b", "c")
    rt <- as.matrix(xtabs(c ~ ., dum))
    dp <- numeric(N)
    t2 <- 2^nt - 1
    nbs <- nb + ns
    nsf <- c()
    for (j1 in 1:nt){
      for (j2 in 1:j1){
        nsf <- rbind(nsf, t(c(j1, j2)))
      }
    }
    colf <- list()

    for (j in 1:(t2)){
      mv <- Dec2Bin(j-1, nt)
      aprt <- as.matrix(apply(rt==0, 1, '==', mv))
      if (t2==1) {
        aprt <- t(aprt)
      }
      dp[colSums(aprt) == nt] <- j
      colf[[j]] <- (1:nt)[mv == 0]
    }
    ndp <- matrix(0, ngs, t2)
    for (i in 1:ngs){
      for (j in 1:t2) {
        ndp[i,j] <- sum(dp == j & grp0 == i)
      }
    }
    Hb <- matrix(0, nb, nb)
    Hbs <- c()
    Jb <- matrix(0, nb, nb)
    Jbs <- c()
    Hs1 <- list()
    Hbs1 <- list()
    Js1 <- list()
    Jbs1 <- list()
    for (i in 1:ngs){
      Hs1[[i]] <- matrix(0,nsg,nsg)
      Hbs1[[i]] <- matrix(0,nb,nsg)
      Js1[[i]] <- matrix(0,nsg,nsg)
      Jbs1[[i]] <- matrix(0,nb,nsg)
      for  (j in 1:t2){
        if (ndp[i,j] != 0){
          if (hetero){
            rj <- rt[dp == j & grp0 == i, colf[[j]], drop=F]
          } else {
            rj <- rt[dp == j, colf[[j]], drop=F]
          }
          nj <- nrow(rj)
          S <- V[[i]][colf[[j]], colf[[j]], drop=F]
          nSj <- nrow(S)

          cholS <- chol(S)
          mat0 <- matrix(0, nSj, nSj)
          dS <- list()
          dS2 <- list()
          if (type == 'UN'){
            for (j1 in 1:nsg){
              dS[[j1]] <- 1 * (S == alp[[i]][j1])
              dS2[[j1]] <- list()
              for (j2 in 1 : nsg){
                dS2[[j1]][[j2]] <- mat0
              }
            }
          }
          if (type == 'CS'){
            dS[[1]] <- diag(nSj)
            dS[[2]] <- matrix(1, nSj, nSj)
            for (j1 in 1:2){
              dS2[[j1]] <- list()
              for (j2 in 1:2){
                dS2[[j1]][[j2]] <- mat0
              }
            }
          }
          if (type == 'AR(1)'){
            dS[[1]] <- S / S[1, 1]
            dS0 <- mat0
            dS[[2]] <- mat0
            for (j1 in 1:2){
              dS2[[j1]] <- list()
              for (j2 in 1:2){
                dS2[[j1]][[j2]] <- mat0
              }
            }
            for (j1 in 1:nSj){
              for (j2 in 1:nSj){
                rp <- abs(j1 - j2)
                dS[[2]][j1, j2] <- rp * alp[[i]][1] * alp[[i]][2] ^ (rp - 1)
                dS2[[1]][[2]][j1, j2] <- rp * alp[[i]][2] ^ (rp - 1)
                dS2[[2]][[1]][j1, j2] <- rp * alp[[i]][2] ^ (rp - 1)
                dS2[[2]][[2]][j1, j2] <- rp * (rp - 1) * alp[[i]][1] *
                  alp[[i]][2] ^ (rp - 2)
              }
            }
          }
          iS <- ginv(S)
          sSr <- backsolve(cholS, t(rj), transpose=T)
          sxb <- list()
          xjb <- list()
          cs_xSr <- list()
          for (jb in 1:nb){
            if (hetero){
              xj <- Xb[[jb]][dp == j & grp0 == i, colf[[j]], drop=F]
            } else {
              xj <- Xb[[jb]][dp == j, colf[[j]], drop=F]
            }
            xjb[[jb]] <- xj
            sxb0 <- backsolve(cholS, t(xj), transpose=T)
            sxb[[jb]] <- sxb0
            cs_xSr[[jb]] <- colSums(sxb0 * sSr)
          }
          qAr <- list()
          cs_rAr <- list()
          rAA <- list()
          AA <- list()
          tris <- numeric(nsg)
          for (js in 1:nsg){
            j1 <- nsf[js, 1]
            j2 <- nsf[js, 2]
            if (sum(colf[[j]] == j1) == 1 & sum(colf[[j]] == j2) == 1){
              AA0 <- -iS %*% dS[[js]] %*% iS
              qr0 <- qr(AA0)
              AA[[js]] <- AA0
              rAA[[js]] <- qr.R(qr0)
              qx0 <- t(qr.Q(qr0)) %*% t(rj)
              rx0 <- qr.R(qr0) %*% t(rj)
              tris[js] <- tr(iS %*% dS[[js]])
              qAr[[js]] <- qx0
              cs_rAr[[js]] <- colSums(qx0 * rx0)
            } else {
              qAr[[js]] <- NA
            }
          }
          for (j1 in 1:nb){
            for (j2 in j1:nb){
              if (!is.na(sxb[[j1]][1]) & !is.na(sxb[[j2]][1])){
                Hb[j1, j2] <- Hb[j1, j2] - sum(sxb[[j1]] * sxb[[j2]])
                Jb[j1, j2] <- Jb[j1, j2] + sum(cs_xSr[[j1]] * cs_xSr[[j2]])
              }
            }
          }
          for (j1 in 1:nsg){
            for (j2 in j1:nsg){
              if (!is.na(qAr[[j1]][1]) & !is.na(qAr[[j2]][1])){
                AA2 <- iS %*% (2 * dS[[j1]] %*% iS %*% dS[[j2]] -
                                 dS2[[j1]][[j2]]) %*% iS
                qr0 <- qr(AA2)
                qx0 <- t(qr.Q(qr0)) %*% t(rj)
                rx0 <- qr.R(qr0) %*% t(rj)
                Hs1[[i]][j1, j2] <- Hs1[[i]][j1, j2] -
                  0.5 * nj * tr(AA[[j1]] %*% dS[[j2]] + iS %*% dS2[[j1]][[j2]]) -
                  0.5 * sum(qx0 * rx0)
                Js1[[i]][j1, j2] <- Js1[[i]][j1, j2] +
                  0.25 * nj * tris[j1] * tris[j2] +
                  0.25 * sum(cs_rAr[[j1]]*cs_rAr[[j2]]) +
                  0.25 * sum(tris[j1] * cs_rAr[[j2]]) +
                  0.25 * sum(tris[j2] * cs_rAr[[j1]])
              }
            }
          }
          for (j1 in 1:nb){
            for (j2 in 1:nsg){
              if (!is.na(qAr[[j2]][1])){
                xrA <- rAA[[j2]] %*% t(xjb[[j1]])
                Hbs1[[i]][j1, j2] <- Hbs1[[i]][j1, j2] + sum(xrA * qAr[[j2]])
                Jbs1[[i]][j1, j2] <- Jbs1[[i]][j1, j2] -
                  0.5 * sum(tris[[j2]] * cs_xSr[[j1]]) -
                  0.5 * sum(cs_xSr[[j1]] * cs_rAr[[j2]])
              }
            }
          }
        }
        for (j1 in 1:nb){
          for (j2 in 1:j1){
            Hb[j1, j2] <- Hb[j2, j1]
            Jb[j1, j2] <- Jb[j2, j1]
          }
        }
        for (j1 in 1:nsg){
          for (j2 in 1:j1){
            Hs1[[i]][j1, j2] <- Hs1[[i]][j2, j1]
            Js1[[i]][j1, j2] <- Js1[[i]][j2, j1]
          }
        }
      }
      Hbs <- cbind(Hbs, Hbs1[[i]])
      Jbs <- cbind(Jbs, Jbs1[[i]])
    }
    Hs <- bdiag(Hs1)
    Js <- bdiag(Js1)
    H <- as.matrix(rbind(cbind(Hb, Hbs), cbind(t(Hbs), Hs)))
    J <- as.matrix(rbind(cbind(Jb, Jbs), cbind(t(Jbs), Js)))
    iII0 <- ginv(-Hb)
    iII <- ginv(-H)
    iIIr0 <- iII0 %*% Jb %*% iII0
    iIIr <- iII %*% J %*% iII
    if (robust) {
      Vtheta <- iIIr
    } else {
      Vtheta <- iII
    }
    res <- c()
    resd <- c()
    for (tp0 in 1:nt){
      dbt <- numeric(nt - 1)
      if (tp0 != 1){
        dbt[tp0 - 1] <- 1
      }
      lsm <- numeric(ng)
      SE_mo <- numeric(ng)
      SE_m <- numeric(ng)
      SE_ro <- numeric(ng)
      SE_r <- numeric(ng)
      ell <- matrix(0, length(beta), ng)
      for (i in 1:ng){
        dbg <- numeric(ng-1)
        if (i != 1){
          dbg[i - 1] <- 1
        }
        dbgt <- numeric((ng - 1) * (nt - 1))
        if (i != 1 & tp0 != 1){
          bgti <- (tp0 - 2) * (ng - 1) + i - 1
          dbgt[bgti] <- 1
        }
        covmeanint <- NULL
        if (!is.null(covariate)){
          if (intcov) {
            covmeanint <- numeric((nt - 1) * length(covmean))
            if (tp0 != 1) {
              covmeanint[seq(tp0 - 1, length(covmeanint),
                             by = (nt - 1))] <- covmean
            }
          }
        }
        ell[, i] <- c(1, dbg, dbt, covmean, dbgt, covmeanint)
        Dt <- c(ell[, i], numeric(ns))
        Vmo <- c(t(ell[, i]) %*% iII0 %*% ell[, i])
        Vm <- c(t(Dt) %*% iII %*% Dt)
        Vro <- c(t(ell[, i]) %*% iIIr0 %*% ell[, i])
        Vr <- c(t(Dt) %*% iIIr %*% Dt)
        SE_mo[i] <- sqrt(Vmo)
        SE_m[i] <- sqrt(Vm)
        SE_ro[i] <- sqrt(Vro)
        SE_r[i] <- sqrt(Vr)
        lsm[i] <- t(ell[, i]) %*% beta
      }
      cbn <- choose(ng, 2)
      comb <- matrix(0, cbn, 2)
      count <- 1
      for (ii in 1:ng){
        for (jj in 1:ng){
          if (ii<jj) {
            comb[count, 1] <- ii
            comb[count, 2] <- jj
            count <- count+1
          }
        }
      }

      delta <- numeric(cbn)
      SE_dmo <- numeric(cbn)
      SE_dm <- numeric(cbn)
      SE_dro <- numeric(cbn)
      SE_dr <- numeric(cbn)
      if (type == "UN"){
        rad1 <- Ncc / (Ncc - nt)
        dad <- Ncc - nt
        if(dad <= 1) {
          stop("Number of complete cases is too small to conduct small
                sample adjustment. Use different ssadjust option.")
        }
        if (rad1 > 1.5 ^ 2) {
          warning("Small sample adjustments have increased the SE by more than
                  1.5 times. May be overly conservative. Might should consider
                  using a different ssadjust option.")
        }
      }
      if (type != "UN"){
        rad1 <- n.data / (n.data - length(beta))
        dad <- (N - ng) * (nt - 1) - n.na
      }
      rad2 <- (N / (N - 1)) * ((n.data - 1) / (n.data - nb))
      for (i in 1:cbn){
        delta[i] <- lsm[comb[i, 2]] - lsm[comb[i, 1]]
        elld <- ell[, comb[i, 2]]-ell[, comb[i, 1]]
        Dt2 <- c(elld, numeric(ns))
        Vdmo <- c(t(elld) %*% iII0 %*% elld)
        Vdm <- c(t(Dt2) %*% iII %*% Dt2)
        Vdro <- c(t(elld) %*% iIIr0 %*% elld)
        Vdr <- c(t(Dt2) %*% iIIr %*% Dt2)
        SE_dmo[i] <- sqrt(Vdmo)
        SE_dm[i] <- sqrt(Vdm)
        SE_dro[i] <- sqrt(Vdro)
        SE_dr[i] <- sqrt(Vdr)
      }
      if (od & !robust) {
        SE <- SE_m
        SEd <- SE_dm
      }
      if (od & robust) {
        SE <- SE_r
        SEd <- SE_dr
      }
      if (!od & !robust) {
        SE <- SE_mo
        SEd <- SE_dmo
      }
      if (!od & robust) {
        SE <- SE_ro
        SEd <- SE_dro
      }
      if (ssadjust == 0) {
        df <- Inf
      }
      if (ssadjust == 1) {
        SE <- SE * sqrt(rad1)
        SEd <- SEd * sqrt(rad1)
        df <- dad
      }
      if (ssadjust == 2) {
        SE <- SE * sqrt(rad2)
        SEd <- SEd * sqrt(rad2)
        df <- dad
      }
      tad <- qt(1 - (1 - conf.level) / 2, df)
      lowerCL <- lsm - tad * SE
      upperCL <- lsm + tad * SE
      lowerCLd <- delta - tad * SEd
      upperCLd <- delta + tad * SEd
      t.value <- delta / SEd
      p.value <- 2 * (1 - pt(abs(delta / SEd), df))
      res0 <- data.frame(time = tl[tp0], group = gl, estimate = lsm,
                         stderr = SE, lowerCL = lowerCL, upperCL = upperCL)
      res0d <- data.frame(time = tl[tp0], group1 = gl[comb[, 2]],
                          group2 = gl[comb[, 1]], estimate = delta,
                          stderr = SEd, lowerCL = lowerCLd, upperCL = upperCLd,
                          t.value = t.value, df = df, p.value = p.value)
      res <- rbind(res, res0)
      resd <- rbind(resd, res0d)
    }
    res <- list(marg.mean = res, diff.mean = resd, Vtheta = Vtheta,
                mmrmObject = fit)
    class(res) <- "mmrmod"
  }
  return(res)
}

