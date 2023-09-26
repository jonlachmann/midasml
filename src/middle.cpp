#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/*
//' @export
// [[Rcpp::export]]
Rcpp::List sglfit_middle(
    int nobs,
    Eigen::MatrixXd &x,
    Eigen::VectorXd &r,
    Eigen::VectorXd b,
    double al,
    double gamma,
    Eigen::VectorXd pf,
    double peps,
    double gw,
    double step
    ) {
    while (true) {
        npass++;
        dif = 0.0;
        if (intr == 1) oldbeta[1] <- b[1]
        if (ni > 0) oldbeta[m[1:ni] + 1] <- b[m[1:ni] + 1]
        pln <- pln + 1
        # ! ----------------- GROUP LOOP -------------------- !
        for (k in seq_len(ngroups)) {
          gidxs <- groups[[k]]
          gw <- sqrt(sum(pf_groups[[k]]))
          skip <- 1
          if (pln == 1) {
            skip <- 0
          }
          if (ju[k] == 1) {
            skip <- 0
          }
          if (skip == 0) {
            # ! --- sg-LASSO PROXIMAL MAP --- !
            prox_res <- prox_sgl(nobs, x_groups[[k]], r, b[gidxs + 1], al, gamma, pf_groups[[k]], peps, gw, steps[k])
            r <- prox_res$r
            b[gidxs + 1] <- prox_res$b
            # ! UPDATE REMAINING VARIABLES
            for (g in gidxs) {
              if (abs(b[g + 1]) > 0) {
                if (pln == 1) {
                  ju[k] <- 1
                }
                d <- oldbeta[g + 1] - b[g + 1]
                dif <- max(dif, maj[g] * d^2)
                if (mm[g] == 0) {
                  ni <- ni + 1
                  if (ni > pmax) break
                  mm[g] <- ni
                  m[ni] <- g # ! RECORD ACTIVE VARIABLES
                }
              }
            }
          } # ! ----------> END GROUP UPDATES
        } # ! ----------> END GROUP LOOP
        if (ni > pmax) break
        if (intr == 1) {
          oldb <- oldbeta[1]
          u <- 0
          while (TRUE) { # ! BEGIN GRADIENT DESCENT
            d <- sum(r) / nobs
            if (d^2 < eps) break
            b[1] <- b[1] + d
            r <- r - d
          } # ! END GRADIENT DESCENT
          d <- b[1] - oldb
          if (abs(d) > 0) dif <- max(dif, d^2)
        }
        if (dif < eps) break
      } # ! ----------> END MIDDLE LOOP
}
*/