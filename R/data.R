#' test dataset for hfun_lr_cal
#'
#' A dataset format as list containing dat.p;dat.g;dat.t;dat.r;dat.m;dat.ng
#'
#' dat.p :data.frame with dim 336 * 56; 336 samples with 4 traits at 14 time points phenotype.
#'
#' dat.g :data.frame with dim 10 * 336; 10 snp genotype of 336 samples.
#'
#' dat.t :numeric =c(1,3,5,7,9,11,13,16,19,23,28,32,39,47)
#'
#' dat.r :numeric number of traits = 4
#'
#' dat.m :logical wether genotype have missing value; TRUE
#'
#' dat.ng :numeric type number of genotype; 2
"test_data"
