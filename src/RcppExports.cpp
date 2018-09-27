// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// c_dtou
List c_dtou(std::vector<std::string> RString, bool rc);
RcppExport SEXP _dtou_c_dtou(SEXP RStringSEXP, SEXP rcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type RString(RStringSEXP);
    Rcpp::traits::input_parameter< bool >::type rc(rcSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dtou(RString, rc));
    return rcpp_result_gen;
END_RCPP
}
// c_dtouDepthLimit
List c_dtouDepthLimit(std::vector<std::string> RString, bool rc, long depth);
RcppExport SEXP _dtou_c_dtouDepthLimit(SEXP RStringSEXP, SEXP rcSEXP, SEXP depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type RString(RStringSEXP);
    Rcpp::traits::input_parameter< bool >::type rc(rcSEXP);
    Rcpp::traits::input_parameter< long >::type depth(depthSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dtouDepthLimit(RString, rc, depth));
    return rcpp_result_gen;
END_RCPP
}
// c_dtouS2
List c_dtouS2(std::vector<std::string> RString, bool rc);
RcppExport SEXP _dtou_c_dtouS2(SEXP RStringSEXP, SEXP rcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type RString(RStringSEXP);
    Rcpp::traits::input_parameter< bool >::type rc(rcSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dtouS2(RString, rc));
    return rcpp_result_gen;
END_RCPP
}
// c_dtouS2DepthLimit
List c_dtouS2DepthLimit(std::vector<std::string> RString, bool rc, long depth);
RcppExport SEXP _dtou_c_dtouS2DepthLimit(SEXP RStringSEXP, SEXP rcSEXP, SEXP depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type RString(RStringSEXP);
    Rcpp::traits::input_parameter< bool >::type rc(rcSEXP);
    Rcpp::traits::input_parameter< long >::type depth(depthSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dtouS2DepthLimit(RString, rc, depth));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dtou_c_dtou", (DL_FUNC) &_dtou_c_dtou, 2},
    {"_dtou_c_dtouDepthLimit", (DL_FUNC) &_dtou_c_dtouDepthLimit, 3},
    {"_dtou_c_dtouS2", (DL_FUNC) &_dtou_c_dtouS2, 2},
    {"_dtou_c_dtouS2DepthLimit", (DL_FUNC) &_dtou_c_dtouS2DepthLimit, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dtou(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
