// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_hfunmap_RCPPEXPORTS_H_GEN_
#define RCPP_hfunmap_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace hfunmap {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("hfunmap", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("hfunmap", "_hfunmap_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in hfunmap");
            }
        }
    }

    inline NumericVector hfun_g_t(Rcpp::NumericVector par, Rcpp::NumericVector Times) {
        typedef SEXP(*Ptr_hfun_g_t)(SEXP,SEXP);
        static Ptr_hfun_g_t p_hfun_g_t = NULL;
        if (p_hfun_g_t == NULL) {
            validateSignature("NumericVector(*hfun_g_t)(Rcpp::NumericVector,Rcpp::NumericVector)");
            p_hfun_g_t = (Ptr_hfun_g_t)R_GetCCallable("hfunmap", "_hfunmap_hfun_g_t");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hfun_g_t(Shield<SEXP>(Rcpp::wrap(par)), Shield<SEXP>(Rcpp::wrap(Times)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline NumericVector hfun_get_mu(arma::vec par, arma::vec Times) {
        typedef SEXP(*Ptr_hfun_get_mu)(SEXP,SEXP);
        static Ptr_hfun_get_mu p_hfun_get_mu = NULL;
        if (p_hfun_get_mu == NULL) {
            validateSignature("NumericVector(*hfun_get_mu)(arma::vec,arma::vec)");
            p_hfun_get_mu = (Ptr_hfun_get_mu)R_GetCCallable("hfunmap", "_hfunmap_hfun_get_mu");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hfun_get_mu(Shield<SEXP>(Rcpp::wrap(par)), Shield<SEXP>(Rcpp::wrap(Times)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline Rcpp::List hfun_get_invsigma(arma::vec par, arma::vec Times, int ranks) {
        typedef SEXP(*Ptr_hfun_get_invsigma)(SEXP,SEXP,SEXP);
        static Ptr_hfun_get_invsigma p_hfun_get_invsigma = NULL;
        if (p_hfun_get_invsigma == NULL) {
            validateSignature("Rcpp::List(*hfun_get_invsigma)(arma::vec,arma::vec,int)");
            p_hfun_get_invsigma = (Ptr_hfun_get_invsigma)R_GetCCallable("hfunmap", "_hfunmap_hfun_get_invsigma");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hfun_get_invsigma(Shield<SEXP>(Rcpp::wrap(par)), Shield<SEXP>(Rcpp::wrap(Times)), Shield<SEXP>(Rcpp::wrap(ranks)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline double hfun_H0_mle(Rcpp::NumericVector par, arma::mat N_phe, Rcpp::NumericVector pheT, int ranks) {
        typedef SEXP(*Ptr_hfun_H0_mle)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_hfun_H0_mle p_hfun_H0_mle = NULL;
        if (p_hfun_H0_mle == NULL) {
            validateSignature("double(*hfun_H0_mle)(Rcpp::NumericVector,arma::mat,Rcpp::NumericVector,int)");
            p_hfun_H0_mle = (Ptr_hfun_H0_mle)R_GetCCallable("hfunmap", "_hfunmap_hfun_H0_mle");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hfun_H0_mle(Shield<SEXP>(Rcpp::wrap(par)), Shield<SEXP>(Rcpp::wrap(N_phe)), Shield<SEXP>(Rcpp::wrap(pheT)), Shield<SEXP>(Rcpp::wrap(ranks)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double hfun_H1_mle(Rcpp::NumericVector par, Rcpp::List pheno, Rcpp::NumericVector Times, int ranks, int ng) {
        typedef SEXP(*Ptr_hfun_H1_mle)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_hfun_H1_mle p_hfun_H1_mle = NULL;
        if (p_hfun_H1_mle == NULL) {
            validateSignature("double(*hfun_H1_mle)(Rcpp::NumericVector,Rcpp::List,Rcpp::NumericVector,int,int)");
            p_hfun_H1_mle = (Ptr_hfun_H1_mle)R_GetCCallable("hfunmap", "_hfunmap_hfun_H1_mle");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hfun_H1_mle(Shield<SEXP>(Rcpp::wrap(par)), Shield<SEXP>(Rcpp::wrap(pheno)), Shield<SEXP>(Rcpp::wrap(Times)), Shield<SEXP>(Rcpp::wrap(ranks)), Shield<SEXP>(Rcpp::wrap(ng)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

}

#endif // RCPP_hfunmap_RCPPEXPORTS_H_GEN_
