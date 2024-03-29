# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393


logidual <- function(y0, k0, rho, cc, lambda, maxiter, cri, crho0) {
    .Call('_glmmkl_logidual', PACKAGE = 'glmmkl', y0, k0, rho, cc, lambda, maxiter, cri, crho0)
}

lrdual <- function(y0, k0, rho, cc, lambda, maxiter, cri, crho0) {
    .Call('_glmmkl_lrdual', PACKAGE = 'glmmkl', y0, k0, rho, cc, lambda, maxiter, cri, crho0)
}

predictspicy <- function(alpha, b, k0) {
    .Call('_glmmkl_predictspicy', PACKAGE = 'glmmkl', alpha, b, k0)
}

mixkercd <- function(xc, xd) {
    .Call('_glmmkl_mixkercd', PACKAGE = 'glmmkl', xc, xd)
}

mixkertestcd <- function(trc, trd, tec, ted) {
    .Call('_glmmkl_mixkertestcd', PACKAGE = 'glmmkl', trc, trd, tec, ted)
}

mixkerc <- function(xc) {
    .Call('_glmmkl_mixkerc', PACKAGE = 'glmmkl', xc)
}

mixkertestc <- function(trc, tec) {
    .Call('_glmmkl_mixkertestc', PACKAGE = 'glmmkl', trc, tec)
}

mixkerd <- function(xd) {
    .Call('_glmmkl_mixkerd', PACKAGE = 'glmmkl', xd)
}

mixkertestd <- function(trd, ted) {
    .Call('_glmmkl_mixkertestd', PACKAGE = 'glmmkl', trd, ted)
}

cvlr <- function(yy, xx, ccsearch, lamsearch, incd, cvwhich, nf, maxiter, cri, crho0) {
    .Call('_glmmkl_cvlr', PACKAGE = 'glmmkl', yy, xx, ccsearch, lamsearch, incd, cvwhich, nf, maxiter, cri, crho0)
}

cvlogi <- function(yy, xx, ccsearch, lamsearch, incd, cvwhich, nf, maxiter, cri, crho0) {
    .Call('_glmmkl_cvlogi', PACKAGE = 'glmmkl', yy, xx, ccsearch, lamsearch, incd, cvwhich, nf, maxiter, cri, crho0)
}

