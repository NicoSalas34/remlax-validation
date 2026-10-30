# ==============================================================================
# Boundary cases: each structure placed on the boundary of its parameter space,
# or near it. Incomplete rank of a us, zero specific variance of a fa, zero
# variance, ar1 or cor correlation close to 1, zero IGE and IEE, overestimated
# rank (rr, fa), extreme scale of y, term with 4 levels.
#
# The criteria stay those of the suite. A component whose true value lies on the
# bound (zero variance) is listed in `bound`: its Hessian standard error has no
# meaning, and its |z| is reported without deciding anything.
#
# Every MME route is tried: automatic, dense, sparse and hybrid inversion.
# ==============================================================================

M_LIM <- list(m_fit("lbfgsb"), m_fit("ai"), m_mme(), m_mme(inversion = "dense"),
              m_mme(inversion = "sparse"), m_mme(inversion = "hybrid"))

boundary <- function(case, bound = character(0), note = NULL) {
  case$family <- "boundary cases"
  case$bound <- bound
  if (!is.null(note)) case$note <- note
  case
}

# --- us of incomplete rank -------------------------------------------------------
lim_us_rank1 <- function(seed) {
  l <- c(0.7, 0.5, -0.4, 0.6)
  Sr <- matrix(0.25, 4, 4); diag(Sr) <- 1
  boundary(case_mt(seed, "boundary_us_rank1_t4", 4, "us", 0L, l %*% t(l), "us", Sr, engines = M_LIM),
         note = "genetic Sigma us of rank 1, residual us: optimum on the boundary of the cone")
}

lim_us8 <- function(seed) {
  set.seed(seed + 7L)
  L <- matrix(stats::rnorm(16, sd = 0.5), 8)
  Sr <- matrix(0.3, 8, 8); diag(Sr) <- 1
  boundary(case_mt(seed, "boundary_us8_near_rank2", 8, "us", 0L, L %*% t(L) + diag(0.01, 8), "us", Sr,
                engines = M_LIM, n_g = 100),
         note = "genetic us(8) nearly of rank 2 (LL' + 0.01 I), as in the absolute benchmark; 100 genotypes")
}

lim_two_us <- function(seed) {
  set.seed(seed)
  t <- 3; D <- mt_design(t, n_g = 150)
  fam <- factor(sprintf("f%02d", (as.integer(D$d$g) - 1L) %/% 5L + 1L))
  Zf <- lapply(D$tn, function(k) { r <- which(D$d$trait == k)
    Matrix::sparseMatrix(i = r, j = as.integer(fam[r]), x = 1, dims = c(D$n, nlevels(fam)),
                         dimnames = list(NULL, levels(fam))) })
  names(Zf) <- D$tn
  Lg <- matrix(c(0.6, 0.4, -0.3, 0, 0.4, 0.3), 3); Sg <- Lg %*% t(Lg)
  lf <- c(0.5, -0.4, 0.3); Sf <- lf %*% t(lf)
  Sr <- matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3)
  V_of <- function(par) V_term(D$Zl, par$sig$g) + V_term(Zf, par$sig$fam) +
    R_section(par$res[[1]], D$tr, D$Cu)
  y <- simulate_y(V_of(list(sig = list(g = Sg, fam = Sf), res = list(Sr))), D$X, seq_len(t))
  m <- rx_model(y, D$X, list(rx_term("g", D$Zl, struct = "us"), rx_term("fam", Zf, struct = "us")),
                rx_residual("us", trait = D$d$trait, unit = D$d$plant))
  truth <- c(sigma_names("g", Sg), sigma_names("fam", Sf), sigma_names("residual", Sr))
  boundary(case_def("boundary_two_us_at_bound", "boundary cases", list(main = list(model = m, V_of = V_of)),
                 M_LIM, truth),
         note = "two us terms on the bound at once: genotype us(3) of rank 2, family us(3) of rank 1")
}

lim_us_cor_strong <- function(seed) {
  R <- matrix(0.97, 3, 3); diag(R) <- 1; sd <- c(0.8, 0.6, 0.7)
  boundary(case_mt(seed, "boundary_us_correlation_097", 3, "us", 0L, diag(sd) %*% R %*% diag(sd),
                "us", matrix(c(1, 0.3, 0.2, 0.3, 0.8, -0.1, 0.2, -0.1, 0.6), 3), engines = M_LIM),
         note = "us(3) of full rank but genetic correlations of 0.97: ill-conditioned interior")
}

# --- fa, rr: zero specific variance, overestimated rank ----------------------------
lim_fa1_psi0 <- function(seed)
  boundary(case_mt(seed, "boundary_fa1_psi_zero", 4, "fa", 1L,
                fa_sigma(matrix(c(0.7, 0.5, -0.4, 0.6)), c(0.2, 0, 0.25, 0.15)),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5)), engines = M_LIM),
         note = "fa rank 1, zero specific variance for trait 2 (Heywood case)")

lim_fa2_rank1 <- function(seed)
  boundary(case_mt(seed, "boundary_fa2_true_rank1", 5, "fa", 2L,
                fa_sigma(matrix(c(0.7, 0.5, -0.3, 0.6, 0.4)), c(0.15, 0.2, 0.2, 0.15, 0.25)),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5, 0.6)), engines = M_LIM, n_g = 100),
         note = "fa rank 2 fitted, truth of rank 1: zero second factor, flat direction")

lim_rr2_rank1 <- function(seed) {
  l <- c(0.7, 0.5, -0.3, 0.6)
  boundary(case_mt(seed, "boundary_rr2_true_rank1", 4, "rr", 2L, l %*% t(l),
                "diag", diag(c(0.8, 0.6, 0.7, 0.5)), engines = M_LIM),
         note = "rr rank 2 fitted, truth of rank 1 without specific variance")
}

# --- zero variances -------------------------------------------------------------------
lim_diag0 <- function(seed)
  boundary(case_mt(seed, "boundary_diag_zero_variance", 3, "diag", 0L, diag(c(0.5, 0, 0.3)),
                "diag", diag(c(1, 0.7, 0.5)), engines = M_LIM),
         bound = "g[2,2]", note = "genetic diag, zero variance for trait 2")

lim_iid0 <- function(seed) {
  set.seed(seed)
  q <- 60; nk <- sample(6:10, q, replace = TRUE)
  g <- factor(rep(sprintf("g%02d", seq_len(q)), nk)); n <- length(g)
  X <- cbind(1, stats::rnorm(n)); Zg <- list(inc(g))
  V_of <- function(par) V_term(Zg, par$sig$g) + R_section(par$res[[1]], NULL, diag(n))
  y <- simulate_y(V_of(list(sig = list(g = 0), res = list(1.0))), X, c(10, 0.5))
  m <- rx_model(y, X, list(rx_term("g", g)), rx_residual())
  boundary(case_def("boundary_iid_zero_variance", "boundary cases", list(main = list(model = m, V_of = V_of)),
                 M_LIM, c(g = 0, residual = 1)),
         bound = "g", note = "one random factor of zero true variance")
}

# --- level correlations close to the bound ---------------------------------------
lim_ar1 <- function(seed, rho, name) {
  set.seed(seed)
  q <- 100; plot <- factor(rep(sprintf("x%03d", 1:q), each = 4)); n <- length(plot)
  X <- matrix(1, n, 1); Zp <- list(inc(plot))
  V_of <- function(par) V_term(Zp, par$sig$plot, C_ar1(q, par$lvl[["plot!ar1"]])) +
    par$res[[1]][1, 1] * diag(n)
  truth <- c(plot = 1, "plot!ar1" = rho, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(plot = 1), res = list(matrix(0.5)), lvl = truth["plot!ar1"])), X, 0)
  m <- rx_model(y, X, list(rx_term("plot", plot, level = "ar1")), rx_residual())
  boundary(case_def(name, "boundary cases", list(main = list(model = m, V_of = V_of)), M_LIM, truth),
         note = sprintf("ar1 term on a transect of 100 plots, rho = %g", rho))
}
lim_ar1_strong <- function(seed) lim_ar1(seed, 0.97, "boundary_ar1_rho_097")
lim_ar1_neg <- function(seed) lim_ar1(seed, -0.6, "boundary_ar1_rho_negative")

lim_res_ar1ar1 <- function(seed) {
  set.seed(seed)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$g <- factor(sample(rep(sprintf("g%03d", 1:100), 6)))
  d$row <- factor(d$row); d$col <- factor(d$col)
  n <- nrow(d); X <- matrix(1, n, 1); Zg <- list(inc(d$g))
  rr <- as.integer(as.character(d$row)); cc <- as.integer(as.character(d$col))
  C_obs <- function(r1, r2) r1^abs(outer(rr, rr, "-")) * r2^abs(outer(cc, cc, "-"))
  V_of <- function(par) V_term(Zg, par$sig$g) +
    R_section(par$res[[1]], NULL, C_obs(par$lvl[["residual!ar1ar11"]], par$lvl[["residual!ar1ar12"]]))
  truth <- c(g = 0.4, residual = 1, "residual!ar1ar11" = 0.95, "residual!ar1ar12" = 0.9)
  y <- simulate_y(V_of(list(sig = list(g = 0.4), res = list(1),
                           lvl = truth[c("residual!ar1ar11", "residual!ar1ar12")])), X, 1)
  d$y <- y
  m <- rx_model(y, X, list(rx_term("g", d$g)), parse_res(~ ar1(row):ar1(col), d))
  boundary(case_def("boundary_residual_ar1xar1_strong", "boundary cases",
                 list(main = list(model = m, V_of = V_of)), M_LIM, truth),
         note = "residual ar1(row):ar1(col) with rho = 0.95 and 0.90 on a 30 x 20 field")
}

lim_cor_strong <- function(seed) {
  set.seed(seed)
  D <- ge_design()
  V_of <- function(par) V_term(D$Zl, par$sig$ge[1, 1] * diag(60), C_cor(4, par$lvl[["ge!cor"]])) +
    par$res[[1]][1, 1] * diag(D$n)
  truth <- c("ge[1,1]" = 0.8, "ge!cor" = 0.97, residual = 0.5)
  y <- simulate_y(V_of(list(sig = list(ge = matrix(0.8)), res = list(matrix(0.5)), lvl = truth["ge!cor"])),
                 D$X, c(1, 0.5, -0.5, 0.2))
  m <- rx_model(y, D$X, list(suppressMessages(rx_term("ge", D$Zl, struct = "iid", level = "cor"))), rx_residual())
  boundary(case_def("boundary_cor_097", "boundary cases", list(main = list(model = m, V_of = V_of)), M_LIM, truth),
         note = "uniform correlation of 0.97 between 4 environments")
}

# --- IGE ---------------------------------------------------------------------------
lim_ige0 <- function(seed) {
  set.seed(seed)
  d <- expand.grid(col = 1:20, row = 1:30)
  d$id <- sprintf("pl%03d", seq_len(nrow(d)))
  d$geno <- sample(rep(sprintf("g%03d", 1:120), 5))
  n <- nrow(d)
  nb <- rx_neighbourhood(coord = cbind(d$row, d$col), group = rep("A", n), level = d$geno, id = d$id,
                         rank = 1, kernel = "none", normalise = TRUE)
  lev <- sort(unique(d$geno))
  Zd <- Matrix::sparseMatrix(i = seq_len(n), j = match(d$geno, lev), x = 1, dims = c(n, 120),
                             dimnames = list(NULL, lev))
  Zn <- nb$level[["A<-A"]][, lev]; dimnames(Zn) <- list(NULL, lev)
  s2 <- c(dge = 0.3, ige = 0, residual = 0.6)
  X <- matrix(1, n, 1)
  V_of <- function(par) V_term(list(Zd), par$sig$dge) + V_term(list(Zn), par$sig$ige) +
    par$res[[1]][1, 1] * diag(n)
  y <- simulate_y(V_of(list(sig = list(dge = 0.3, ige = 0), res = list(matrix(0.6)))), X, 10)
  m <- rx_model(y, X, list(rx_term("dge", list(Zd)), rx_term("ige", list(Zn))), rx_residual())
  boundary(case_def("boundary_ige_zero", "boundary cases", list(main = list(model = m, V_of = V_of)), M_LIM, s2),
         bound = "ige", note = "monoculture 30 x 20, DGE iid + IGE iid of zero true variance")
}

lim_mixture <- function(seed) {
  a <- c(1, 0); b <- c(-0.6, 0.8); cc <- c(0.6, 0.8)
  R <- rbind(a, b, cc) %*% t(rbind(a, b, cc))
  boundary(case_ige_mixture(seed, name = "boundary_mixture_iee_zero_us_rank2", share_iee = 0, Rg = R),
         bound = "iee",
         note = paste("mixture of two species: IEE of zero variance, and (DGE, IGE intra, IGE inter)",
                      "us(3) of rank 2 for each species"))
}

# --- scale of y and few levels ----------------------------------------------------
lim_scale <- function(seed, f, name, n_env = 12) {
  set.seed(seed)
  d <- expand.grid(g = sprintf("g%02d", 1:50), env = sprintf("e%02d", seq_len(n_env)))
  d <- d[sort(sample(nrow(d), round(0.92 * nrow(d)))), ]
  n <- nrow(d); X <- matrix(1, n, 1)
  Zg <- list(inc(d$g)); Ze <- list(inc(d$env))
  V_of <- function(par) V_term(Zg, par$sig$g) + V_term(Ze, par$sig$env) + par$res[[1]][1, 1] * diag(n)
  truth <- c(g = 0.5, env = 0.8, residual = 1.0) * f^2
  y <- simulate_y(V_of(list(sig = list(g = truth[["g"]], env = truth[["env"]]), res = list(matrix(truth[["residual"]])))),
                 X, 5 * f)
  m <- rx_model(y, X, list(rx_term("g", d$g), rx_term("env", d$env)), rx_residual())
  boundary(case_def(name, "boundary cases", list(main = list(model = m, V_of = V_of)), M_LIM, truth),
         note = sprintf("genotype x environment crossed, %d environments, y on scale %g", n_env, f))
}
lim_scale_large <- function(seed) lim_scale(seed, 1e3, "boundary_scale_1e3")
lim_scale_small <- function(seed) lim_scale(seed, 1e-3, "boundary_scale_1em3")
lim_few_levels <- function(seed) lim_scale(seed, 1, "boundary_term_4_levels", n_env = 4)

CATALOGUE_BOUNDARY <- list(
  boundary_us_rank1_t4 = lim_us_rank1, boundary_us8_near_rank2 = lim_us8,
  boundary_two_us_at_bound = lim_two_us, boundary_us_correlation_097 = lim_us_cor_strong,
  boundary_fa1_psi_zero = lim_fa1_psi0, boundary_fa2_true_rank1 = lim_fa2_rank1,
  boundary_rr2_true_rank1 = lim_rr2_rank1, boundary_diag_zero_variance = lim_diag0,
  boundary_iid_zero_variance = lim_iid0, boundary_ar1_rho_097 = lim_ar1_strong,
  boundary_ar1_rho_negative = lim_ar1_neg, boundary_residual_ar1xar1_strong = lim_res_ar1ar1,
  boundary_cor_097 = lim_cor_strong, boundary_ige_zero = lim_ige0,
  boundary_mixture_iee_zero_us_rank2 = lim_mixture, boundary_scale_1e3 = lim_scale_large,
  boundary_scale_1em3 = lim_scale_small, boundary_term_4_levels = lim_few_levels)
CATALOGUE <- c(CATALOGUE, CATALOGUE_BOUNDARY)
