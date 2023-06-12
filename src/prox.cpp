#include <RcppEigen.h>
#include <queue>
#include <random>
// [[Rcpp::depends(RcppEigen)]]

//' @export
// [[Rcpp::export]]
Rcpp::List prox_sgl(
    int nobs,
    Eigen::MatrixXd x,
    Eigen::VectorXd r,
    Eigen::VectorXd b,
    double al,
    double gamma,
    Eigen::VectorXd pf,
    double peps,
    double gw,
    double step
    ) {

    Eigen::VectorXd bold = b;
    double big = 9.9E30;
    double vg_const = step * gw * al * (1 - gamma);
    //#! -------- BEGIN PROGRAM -------- #!
    while (true) {
        bold = b;
        //#!--------- LASSO PART ----------#!
        Eigen::VectorXd u = b + step * x.transpose() * r / nobs;
        Eigen::VectorXd v = u.cwiseAbs() - step * al * gamma * pf;
        v = v.cwiseMax(Eigen::VectorXd::Zero(v.size()));
        b = v.array() * u.array().sign();

        //#!--------- g-LASSO PART ----------#!
        //#! L2 norm of b_g
        double normg = sqrt(b.transpose() * b);
        double vg = big;
        if (normg != 0.0) {
            vg = vg_const / normg;
        }

        Eigen::VectorXd scl = (1.0 - pf.array() * vg).cwiseMax(0.0);
        b = b * scl;
        Eigen::VectorXd d = b - bold;
        r = r - x * d;
        double maxg = d.lpNorm<Eigen::Infinity>();

        //#!--------- CHECK CONVERGENCE ----------#!
        if (maxg < peps) break;
    }
    Rcpp::List result;
    result["b"] = b;
    result["r"] = r;
    return result;
}