#######################################################################
# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' uirtestimate
#' @param data -
#' @param model_ -
#' @return list
#' @keywords internal
uirtestimate <- function(data, model_,convergenceEpsilon) {
    .Call('IRTpp_uirtestimate', PACKAGE = 'IRTpp', data, model_,convergenceEpsilon)
}

#######################################################################
#'  @name abilityinterface
#'  @param zita_par 
#'  @param data 
#'  @param model_ 
#'  @param method 
#'  @param matrix_flag 
#'  @param prob_matrix 
#' @return list
#' @keywords internal
abilityinterface <- function(zita_par, data, model_, method, matrix_flag, prob_matrix) {
    .Call('IRTpp_abilityinterface', PACKAGE = 'IRTpp', zita_par, data, model_, method, matrix_flag, prob_matrix)
}

#######################################################################
#' mapinterface
#'  @param zita_par -
#'  @param dat -
#'  @param e_model -
#'  @param matrix_flag -
#'  @param prob_matrix -
#' @return list
#' @keywords internal
mapinterface <- function(zita_par, dat, e_model, matrix_flag, prob_matrix) {
    .Call('IRTpp_mapinterface', PACKAGE = 'IRTpp', zita_par, dat, e_model, matrix_flag, prob_matrix)
}

#######################################################################
#' eapinterface
#'  @param zita_par  -
#'  @param dat -
#'  @param e_model -
#'  @param matrix_flag -
#'  @param prob_matrix -
#' @return list
#' @keywords internal
eapinterface <- function(zita_par, dat, e_model, matrix_flag, prob_matrix) {
    .Call('IRTpp_eapinterface', PACKAGE = 'IRTpp', zita_par, dat, e_model, matrix_flag, prob_matrix)
}

