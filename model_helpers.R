# Model definitions and diagnostic probability functions.

# ------------------------------------------------------------
# TWO-TEST, THREE-POPULATION MODEL
# test 1 = serology; test 2 = coprology
# ------------------------------------------------------------

ferik_model_string <- function(dependence = TRUE) {
  if (dependence) {
    return("
model {
  for (i in 1:3) {
    y[i,1:4] ~ dmulti(pr[i,1:4], n[i])

    # observed order = --, +-, -+, ++
    pr[i,1] <- p[i] * d00 + (1-p[i]) * nd00
    pr[i,2] <- p[i] * d10 + (1-p[i]) * nd10
    pr[i,3] <- p[i] * d01 + (1-p[i]) * nd01
    pr[i,4] <- p[i] * d11 + (1-p[i]) * nd11
  }

  # Disease-positive cell probabilities
  d11 <- SeS * SeC + covSe
  d10 <- SeS * (1-SeC) - covSe
  d01 <- (1-SeS) * SeC - covSe
  d00 <- (1-SeS) * (1-SeC) + covSe

  # Disease-negative cell probabilities
  nd11 <- (1-SpS) * (1-SpC) + covSp
  nd10 <- (1-SpS) * SpC - covSp
  nd01 <- SpS * (1-SpC) - covSp
  nd00 <- SpS * SpC + covSp

  # Frechet bounds for covariance
  loSe <- max(-SeS*SeC, -(1-SeS)*(1-SeC))
  hiSe <- min(SeS*(1-SeC), (1-SeS)*SeC)
  uSe ~ dunif(0,1)
  covSe <- loSe + uSe * (hiSe-loSe)

  loSp <- max(-SpS*SpC, -(1-SpS)*(1-SpC))
  hiSp <- min(SpS*(1-SpC), (1-SpS)*SpC)
  uSp ~ dunif(0,1)
  covSp <- loSp + uSp * (hiSp-loSp)

  # Priors
  p[1] ~ dbeta(ap1, bp1)
  p[2] ~ dbeta(ap2, bp2)
  p[3] ~ dbeta(ap3, bp3)

  SeS ~ dbeta(aSeS, bSeS)
  SpS ~ dbeta(aSpS, bSpS)
  SeC ~ dbeta(aSeC, bSeC)
  SpC ~ dbeta(aSpC, bSpC)
}
")
  }

  "
model {
  for (i in 1:3) {
    y[i,1:4] ~ dmulti(pr[i,1:4], n[i])

    pr[i,1] <- p[i] * d00 + (1-p[i]) * nd00
    pr[i,2] <- p[i] * d10 + (1-p[i]) * nd10
    pr[i,3] <- p[i] * d01 + (1-p[i]) * nd01
    pr[i,4] <- p[i] * d11 + (1-p[i]) * nd11
  }

  d11 <- SeS * SeC
  d10 <- SeS * (1-SeC)
  d01 <- (1-SeS) * SeC
  d00 <- (1-SeS) * (1-SeC)

  nd11 <- (1-SpS) * (1-SpC)
  nd10 <- (1-SpS) * SpC
  nd01 <- SpS * (1-SpC)
  nd00 <- SpS * SpC

  p[1] ~ dbeta(ap1, bp1)
  p[2] ~ dbeta(ap2, bp2)
  p[3] ~ dbeta(ap3, bp3)

  SeS ~ dbeta(aSeS, bSeS)
  SpS ~ dbeta(aSpS, bSpS)
  SeC ~ dbeta(aSeC, bSeC)
  SpC ~ dbeta(aSpC, bSpC)
}
"
}

make_ferik_data <- function(prior, dependence = TRUE) {
  list(
    y = ferik_counts,
    n = as.integer(ferik_n),
    ap1 = prior$p_arab[1], bp1 = prior$p_arab[2],
    ap2 = prior$p_goran[1], bp2 = prior$p_goran[2],
    ap3 = prior$p_fulani[1], bp3 = prior$p_fulani[2],
    aSeS = prior$se_sero[1], bSeS = prior$se_sero[2],
    aSpS = prior$sp_sero[1], bSpS = prior$sp_sero[2],
    aSeC = prior$se_copro[1], bSeC = prior$se_copro[2],
    aSpC = prior$sp_copro[1], bSpC = prior$sp_copro[2]
  )
}

ferik_monitor <- function(dependence = TRUE) {
  base <- c("p", "SeS", "SpS", "SeC", "SpC")
  if (dependence) c(base, "covSe", "covSp") else base
}

# ------------------------------------------------------------
# THREE-TEST, ONE-POPULATION MODEL
# Explicit names:
# L = liver inspection, C = coprology, S = serology
#
# Pairwise covariance expansion with third-order interaction = 0.
# Invalid joint-probability states are rejected using Bernoulli
# validity constraints.
# ------------------------------------------------------------

slaughter_model_string <- function(dependence = TRUE, constrain_sp_liver = FALSE) {

  sp_liver_prior <- if (constrain_sp_liver) {
    "SpL ~ dbeta(aSpL, bSpL) T(0.9000001, 1)"
  } else {
    "SpL ~ dbeta(aSpL, bSpL)"
  }

  if (!dependence) {
    return(sprintf("
model {
  y[1:8] ~ dmulti(pr[1:8], n)

  # Disease-positive cells, order L/C/S:
  # ---, +--, -+-, ++-, --+, +-+, -++, +++
  d[1] <- (1-SeL)*(1-SeC)*(1-SeS)
  d[2] <- SeL*(1-SeC)*(1-SeS)
  d[3] <- (1-SeL)*SeC*(1-SeS)
  d[4] <- SeL*SeC*(1-SeS)
  d[5] <- (1-SeL)*(1-SeC)*SeS
  d[6] <- SeL*(1-SeC)*SeS
  d[7] <- (1-SeL)*SeC*SeS
  d[8] <- SeL*SeC*SeS

  # Disease-negative cells use false-positive probabilities
  fL <- 1-SpL
  fC <- 1-SpC
  fS <- 1-SpS

  nd[1] <- (1-fL)*(1-fC)*(1-fS)
  nd[2] <- fL*(1-fC)*(1-fS)
  nd[3] <- (1-fL)*fC*(1-fS)
  nd[4] <- fL*fC*(1-fS)
  nd[5] <- (1-fL)*(1-fC)*fS
  nd[6] <- fL*(1-fC)*fS
  nd[7] <- (1-fL)*fC*fS
  nd[8] <- fL*fC*fS

  for (k in 1:8) {
    pr[k] <- p*d[k] + (1-p)*nd[k]
  }

  p ~ dbeta(ap, bp)
  SeL ~ dbeta(aSeL, bSeL)
  %s
  SeC ~ dbeta(aSeC, bSeC)
  SpC ~ dbeta(aSpC, bSpC)
  SeS ~ dbeta(aSeS, bSeS)
  SpS ~ dbeta(aSpS, bSpS)
}
", sp_liver_prior))
  }

  sprintf("
model {
  y[1:8] ~ dmulti(pr[1:8], n)

  eps <- 1.0E-10

  # --------------------------
  # Disease-positive component
  # --------------------------
  mL <- SeL
  mC <- SeC
  mS <- SeS

  loSeLC <- max(-mL*mC, -(1-mL)*(1-mC))
  hiSeLC <- min(mL*(1-mC), (1-mL)*mC)
  uSeLC ~ dunif(0,1)
  covSeLC <- loSeLC + uSeLC*(hiSeLC-loSeLC)

  loSeLS <- max(-mL*mS, -(1-mL)*(1-mS))
  hiSeLS <- min(mL*(1-mS), (1-mL)*mS)
  uSeLS ~ dunif(0,1)
  covSeLS <- loSeLS + uSeLS*(hiSeLS-loSeLS)

  loSeCS <- max(-mC*mS, -(1-mC)*(1-mS))
  hiSeCS <- min(mC*(1-mS), (1-mC)*mS)
  uSeCS ~ dunif(0,1)
  covSeCS <- loSeCS + uSeCS*(hiSeCS-loSeCS)

  eLC_D <- mL*mC + covSeLC
  eLS_D <- mL*mS + covSeLS
  eCS_D <- mC*mS + covSeCS
  eLCS_D <- mL*mC*mS + covSeLC*mS + covSeLS*mC + covSeCS*mL

  dRaw[8] <- eLCS_D
  dRaw[4] <- eLC_D - eLCS_D
  dRaw[6] <- eLS_D - eLCS_D
  dRaw[7] <- eCS_D - eLCS_D
  dRaw[2] <- mL - eLC_D - eLS_D + eLCS_D
  dRaw[3] <- mC - eLC_D - eCS_D + eLCS_D
  dRaw[5] <- mS - eLS_D - eCS_D + eLCS_D
  dRaw[1] <- 1 - mL - mC - mS + eLC_D + eLS_D + eCS_D - eLCS_D

  validD <- step(dRaw[1]-eps) * step(dRaw[2]-eps) * step(dRaw[3]-eps) *
            step(dRaw[4]-eps) * step(dRaw[5]-eps) * step(dRaw[6]-eps) *
            step(dRaw[7]-eps) * step(dRaw[8]-eps)
  oneD ~ dbern(validD)

  for (k in 1:8) {
    dPos[k] <- max(eps, dRaw[k])
  }
  sumD <- sum(dPos[])
  for (k in 1:8) {
    d[k] <- dPos[k] / sumD
  }

  # --------------------------
  # Disease-negative component
  # false-positive marginals; covariance is unchanged by complementing
  # both indicators, so manuscript Sp covariance terms apply here.
  # --------------------------
  fL <- 1-SpL
  fC <- 1-SpC
  fS <- 1-SpS

  loSpLC <- max(-fL*fC, -(1-fL)*(1-fC))
  hiSpLC <- min(fL*(1-fC), (1-fL)*fC)
  uSpLC ~ dunif(0,1)
  covSpLC <- loSpLC + uSpLC*(hiSpLC-loSpLC)

  loSpLS <- max(-fL*fS, -(1-fL)*(1-fS))
  hiSpLS <- min(fL*(1-fS), (1-fL)*fS)
  uSpLS ~ dunif(0,1)
  covSpLS <- loSpLS + uSpLS*(hiSpLS-loSpLS)

  loSpCS <- max(-fC*fS, -(1-fC)*(1-fS))
  hiSpCS <- min(fC*(1-fS), (1-fC)*fS)
  uSpCS ~ dunif(0,1)
  covSpCS <- loSpCS + uSpCS*(hiSpCS-loSpCS)

  eLC_N <- fL*fC + covSpLC
  eLS_N <- fL*fS + covSpLS
  eCS_N <- fC*fS + covSpCS
  eLCS_N <- fL*fC*fS + covSpLC*fS + covSpLS*fC + covSpCS*fL

  nRaw[8] <- eLCS_N
  nRaw[4] <- eLC_N - eLCS_N
  nRaw[6] <- eLS_N - eLCS_N
  nRaw[7] <- eCS_N - eLCS_N
  nRaw[2] <- fL - eLC_N - eLS_N + eLCS_N
  nRaw[3] <- fC - eLC_N - eCS_N + eLCS_N
  nRaw[5] <- fS - eLS_N - eCS_N + eLCS_N
  nRaw[1] <- 1 - fL - fC - fS + eLC_N + eLS_N + eCS_N - eLCS_N

  validN <- step(nRaw[1]-eps) * step(nRaw[2]-eps) * step(nRaw[3]-eps) *
            step(nRaw[4]-eps) * step(nRaw[5]-eps) * step(nRaw[6]-eps) *
            step(nRaw[7]-eps) * step(nRaw[8]-eps)
  oneN ~ dbern(validN)

  for (k in 1:8) {
    nPos[k] <- max(eps, nRaw[k])
  }
  sumN <- sum(nPos[])
  for (k in 1:8) {
    nd[k] <- nPos[k] / sumN
  }

  for (k in 1:8) {
    pr[k] <- p*d[k] + (1-p)*nd[k]
  }

  # Priors
  p ~ dbeta(ap, bp)
  SeL ~ dbeta(aSeL, bSeL)
  %s
  SeC ~ dbeta(aSeC, bSeC)
  SpC ~ dbeta(aSpC, bSpC)
  SeS ~ dbeta(aSeS, bSeS)
  SpS ~ dbeta(aSpS, bSpS)
}
", sp_liver_prior)
}

make_slaughter_data <- function(prior) {
  list(
    y = as.integer(slaughter_counts),
    n = as.integer(sum(slaughter_counts)),
    oneD = 1L,
    oneN = 1L,
    ap = prior$p[1], bp = prior$p[2],
    aSeL = prior$se_liver[1], bSeL = prior$se_liver[2],
    aSpL = prior$sp_liver[1], bSpL = prior$sp_liver[2],
    aSeC = prior$se_copro[1], bSeC = prior$se_copro[2],
    aSpC = prior$sp_copro[1], bSpC = prior$sp_copro[2],
    aSeS = prior$se_sero[1], bSeS = prior$se_sero[2],
    aSpS = prior$sp_sero[1], bSpS = prior$sp_sero[2]
  )
}

slaughter_monitor <- function(dependence = TRUE) {
  base <- c("p", "SeL", "SpL", "SeC", "SpC", "SeS", "SpS")
  if (dependence) {
    c(base, "covSeLC", "covSeLS", "covSeCS", "covSpLC", "covSpLS", "covSpCS")
  } else {
    base
  }
}


