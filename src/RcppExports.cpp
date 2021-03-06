// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/hfunmap.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// hfun_g_t
NumericVector hfun_g_t(Rcpp::NumericVector par, Rcpp::NumericVector Times);
static SEXP _hfunmap_hfun_g_t_try(SEXP parSEXP, SEXP TimesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Times(TimesSEXP);
    rcpp_result_gen = Rcpp::wrap(hfun_g_t(par, Times));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _hfunmap_hfun_g_t(SEXP parSEXP, SEXP TimesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_hfunmap_hfun_g_t_try(parSEXP, TimesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hfun_get_mu
NumericVector hfun_get_mu(arma::vec par, arma::vec Times);
static SEXP _hfunmap_hfun_get_mu_try(SEXP parSEXP, SEXP TimesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Times(TimesSEXP);
    rcpp_result_gen = Rcpp::wrap(hfun_get_mu(par, Times));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _hfunmap_hfun_get_mu(SEXP parSEXP, SEXP TimesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_hfunmap_hfun_get_mu_try(parSEXP, TimesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hfun_get_invsigma
Rcpp::List hfun_get_invsigma(arma::vec par, arma::vec Times, int ranks);
static SEXP _hfunmap_hfun_get_invsigma_try(SEXP parSEXP, SEXP TimesSEXP, SEXP ranksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Times(TimesSEXP);
    Rcpp::traits::input_parameter< int >::type ranks(ranksSEXP);
    rcpp_result_gen = Rcpp::wrap(hfun_get_invsigma(par, Times, ranks));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _hfunmap_hfun_get_invsigma(SEXP parSEXP, SEXP TimesSEXP, SEXP ranksSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_hfunmap_hfun_get_invsigma_try(parSEXP, TimesSEXP, ranksSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hfun_H0_mle
double hfun_H0_mle(Rcpp::NumericVector par, arma::mat N_phe, Rcpp::NumericVector pheT, int ranks);
static SEXP _hfunmap_hfun_H0_mle_try(SEXP parSEXP, SEXP N_pheSEXP, SEXP pheTSEXP, SEXP ranksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type N_phe(N_pheSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pheT(pheTSEXP);
    Rcpp::traits::input_parameter< int >::type ranks(ranksSEXP);
    rcpp_result_gen = Rcpp::wrap(hfun_H0_mle(par, N_phe, pheT, ranks));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _hfunmap_hfun_H0_mle(SEXP parSEXP, SEXP N_pheSEXP, SEXP pheTSEXP, SEXP ranksSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_hfunmap_hfun_H0_mle_try(parSEXP, N_pheSEXP, pheTSEXP, ranksSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hfun_H1_mle
double hfun_H1_mle(Rcpp::NumericVector par, Rcpp::List pheno, Rcpp::NumericVector Times, int ranks, int ng);
static SEXP _hfunmap_hfun_H1_mle_try(SEXP parSEXP, SEXP phenoSEXP, SEXP TimesSEXP, SEXP ranksSEXP, SEXP ngSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pheno(phenoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Times(TimesSEXP);
    Rcpp::traits::input_parameter< int >::type ranks(ranksSEXP);
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    rcpp_result_gen = Rcpp::wrap(hfun_H1_mle(par, pheno, Times, ranks, ng));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _hfunmap_hfun_H1_mle(SEXP parSEXP, SEXP phenoSEXP, SEXP TimesSEXP, SEXP ranksSEXP, SEXP ngSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_hfunmap_hfun_H1_mle_try(parSEXP, phenoSEXP, TimesSEXP, ranksSEXP, ngSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _hfunmap_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("NumericVector(*hfun_g_t)(Rcpp::NumericVector,Rcpp::NumericVector)");
        signatures.insert("NumericVector(*hfun_get_mu)(arma::vec,arma::vec)");
        signatures.insert("Rcpp::List(*hfun_get_invsigma)(arma::vec,arma::vec,int)");
        signatures.insert("double(*hfun_H0_mle)(Rcpp::NumericVector,arma::mat,Rcpp::NumericVector,int)");
        signatures.insert("double(*hfun_H1_mle)(Rcpp::NumericVector,Rcpp::List,Rcpp::NumericVector,int,int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _hfunmap_RcppExport_registerCCallable() { 
    R_RegisterCCallable("hfunmap", "_hfunmap_hfun_g_t", (DL_FUNC)_hfunmap_hfun_g_t_try);
    R_RegisterCCallable("hfunmap", "_hfunmap_hfun_get_mu", (DL_FUNC)_hfunmap_hfun_get_mu_try);
    R_RegisterCCallable("hfunmap", "_hfunmap_hfun_get_invsigma", (DL_FUNC)_hfunmap_hfun_get_invsigma_try);
    R_RegisterCCallable("hfunmap", "_hfunmap_hfun_H0_mle", (DL_FUNC)_hfunmap_hfun_H0_mle_try);
    R_RegisterCCallable("hfunmap", "_hfunmap_hfun_H1_mle", (DL_FUNC)_hfunmap_hfun_H1_mle_try);
    R_RegisterCCallable("hfunmap", "_hfunmap_RcppExport_validate", (DL_FUNC)_hfunmap_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_hfunmap_hfun_g_t", (DL_FUNC) &_hfunmap_hfun_g_t, 2},
    {"_hfunmap_hfun_get_mu", (DL_FUNC) &_hfunmap_hfun_get_mu, 2},
    {"_hfunmap_hfun_get_invsigma", (DL_FUNC) &_hfunmap_hfun_get_invsigma, 3},
    {"_hfunmap_hfun_H0_mle", (DL_FUNC) &_hfunmap_hfun_H0_mle, 4},
    {"_hfunmap_hfun_H1_mle", (DL_FUNC) &_hfunmap_hfun_H1_mle, 5},
    {"_hfunmap_RcppExport_registerCCallable", (DL_FUNC) &_hfunmap_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_hfunmap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
