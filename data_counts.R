# ============================================================
# Chad Fasciola BLCM - final aggregate data counts and priors
# ============================================================

# Herd counts
# columns: serology/copro = --, +-, -+, ++
ferik_counts <- rbind(
  Arab   = c(67, 2, 6, 1),
  Goran  = c(22, 0, 1, 0),
  Fulani = c(32, 78, 1, 7)
)
colnames(ferik_counts) <- c("sero-_copro-", "sero+_copro-", "sero-_copro+", "sero+_copro+")
ferik_n <- rowSums(ferik_counts)

stopifnot(identical(as.integer(ferik_n), c(76L, 23L, 118L)))
stopifnot(sum(ferik_counts) == 217)

# Slaughterhouse sheep/goat counts after exclusion of the bovine-coded record.
# order: liver/copro/sero = ---, +--, -+-, ++-, --+, +-+, -++, +++
slaughter_counts <- c(
  86,  # ---
  5,   # +--
  4,   # -+-
  0,   # ++-
  9,   # --+
  1,   # +-+
  3,   # -++
  1    # +++
)
names(slaughter_counts) <- c("---", "+--", "-+-", "++-", "--+", "+-+", "-++", "+++")
stopifnot(sum(slaughter_counts) == 109)

# ----------------------------
# Baseline prior constants
# ----------------------------

ferik_prior <- list(
  p_arab   = c(3.901, 11.037),
  p_goran  = c(1.148, 80.948),
  p_fulani = c(10.071, 11.229),
  se_sero  = c(53.581, 2.626),
  sp_sero  = c(54.706, 5.042),
  se_copro = c(4.842, 3.561),
  sp_copro = c(14.844, 4.461)
)

# In the slaughterhouse model, the ELISA and coprology test priors are replaced
# by moment-matched Beta distributions from the primary herd posterior.
# The following constants specify the slaughterhouse prevalence prior
# and the Beta(1,1) liver-inspection priors.
slaughter_prior <- list(
  p        = c(3.693, 10.016),
  se_liver = c(1, 1),
  sp_liver = c(1, 1),
  se_copro = c(1, 1),
  sp_copro = c(1, 1),
  se_sero  = c(1, 1),
  sp_sero  = c(1, 1)
)

widen_beta <- function(x, factor = 5) {
  stopifnot(length(x) == 2, factor > 0)
  x / factor
}

uninformative_beta <- c(1, 1)

