//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::interfaces(r,cpp)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace arma;


//' Single trait logistic curve fit
//'
//' @param par  numeric of curve par
//' @param Times numeric of times
//'
//' @return vector of logistic curve number
//'
//' @export
//[[Rcpp::export]]
NumericVector hfun_g_t(Rcpp::NumericVector par, Rcpp::NumericVector Times){
  Rcpp::NumericVector u1 = wrap(par[0]/(1 + par[1] * exp(-par[2]*Times)));
  return u1;
}
//' Multi-trait logistic curve fit
//' @param par numeric of curve par
//' @param Times numeric of times
//'
//' @return  matrix dim = 1*length_n.trait*length_n.times
//'
//' @export
//[[Rcpp::export]]
NumericVector hfun_get_mu(arma::vec par, arma::vec Times){
  int rank = par.size()/3;
  int n = Times.size();
  mat u(n,rank);
  for(int i = 0; i < rank; ++i){
    u.col(i) = arma::vec (par[0+3*i]/(1 + par[1+3*i] * exp(-par[2+3*i]*Times)));
  }
  NumericMatrix output = wrap(u.as_row());
  output.attr("dim") = R_NilValue;
  return output;
}
//' Multi-trait SAD1 inv & det
//' @param par numeric of covariance matrix par
//' @param Times numeric of times
//' @param ranks int, number of traits
//' @return  list $sigma_1=sigma_1;$det=det_sigma_1
//'
//' @export
//[[Rcpp::export]]
Rcpp::List hfun_get_invsigma(arma::vec par, arma::vec Times, int ranks){
  int m = Times.size();
  int n = m * ranks;
  arma::mat sigma_1(n,n,fill::zeros);
  double sigma_1_det;
  if(ranks==1){
    double phi = par[0];
    double v2 = par[1];
    arma::vec diag1(m,fill::ones);
    sigma_1.diag(0) = diag1;
    for(int i=0;i<n-1;i++){
      sigma_1(i,i+1) = -pow(phi,(Times[i+1] - Times[i]));
      sigma_1(i,i) = 1 + pow(phi,(2 * (Times[i+1] - Times[i])));
      sigma_1(i+1,i) = sigma_1(i,i+1);
    }
    sigma_1 = sigma_1 * abs(1/v2);
    sigma_1_det = pow(1/v2,m);

  }
  else if(ranks==2){
    double phi1 = par[0];
    double phi2 = par[1];
    double r1 = par[2];
    double r2 = par[3];
    double rho = par[4];
    for(int i = 0;i<m-1;i++){
      sigma_1.submat(span(2*i,2*i+1),span(2*i+2,2*i+3)) = arma::mat {{-pow(phi1,(Times[i+1]-Times[i]))/(pow(r1,2) * (1 - pow(rho,2))),
                                                            rho * pow(phi2,(Times[i+1]-Times[i]))/(r1 * r2 * (1 - pow(rho,2)))},
                                                            {rho * pow(phi1,(Times[i+1]-Times[i]))/(r1 * r2 * (1 - pow(rho,2))),
                                                             -pow(phi2,(Times[i+1]-Times[i]))/(pow(r2,2) * (1 - pow(rho,2)))}};
      sigma_1.submat(span(2*i,2*i+1),span(2*i,2*i+1)) = arma::mat {{1/(pow(r1,2)*(1-pow(rho,2)))*(1+pow(phi1,(Times[i+1]-Times[i]))*pow(phi1,(Times[i+1]-Times[i]))),
                                                          -rho/(r1*r2*(1-pow(rho,2)))*(1+pow(phi2,(Times[i+1]-Times[i]))*pow(phi1,(Times[i+1]-Times[i])))},
                                                          {-rho/(r1*r2*(1-pow(rho,2)))*(1+pow(phi1,(Times[i+1]-Times[i]))*pow(phi2,(Times[i+1]-Times[i]))),
                                                           1/(pow(r2,2)*(1-pow(rho,2)))*(1+pow(phi2,(Times[i+1]-Times[i]))*pow(phi2,(Times[i+1]-Times[i])))}};
      sigma_1.submat(span(2*i+2,2*i+3),span(2*i,2*i+1)) = sigma_1.submat(span(2*i,2*i+1),span(2*i+2,2*i+3)).t();
    }
    sigma_1.submat(span(n-2,n-1),span(n-2,n-1)) = arma::mat {{1/(pow(r1,2)*(1-pow(rho,2))),-rho/(r1*r2*(1-pow(rho,2)))},
                                               {-rho/(r1*r2*(1-pow(rho,2))),1/(pow(r2,2)*(1-pow(rho,2)))}};
    sigma_1_det = pow(1/(pow(r1,2)*pow(r2,2)*(1-pow(rho,2))),m);

  }
  else if(ranks==4){
    double phi1 = par[0];
    double phi2 = par[1];
    double phi3 = par[2];
    double phi4 = par[3];
    double r1 = par[4];
    double r2 = par[5];
    double r3 = par[6];
    double r4 = par[7];
    double ro1 = par[8];
    double ro2 = par[9];
    double ro3 = par[10];
    double ro4 = par[11];
    double ro5 = par[12];
    double ro6 = par[13];
    arma::mat v = { {phi1,0,0,0},{0,phi2,0,0},{0,0,phi3,0},{0,0,0,phi4} };
    arma::mat s = {{r1*r1,r1*r2*ro1,r1*r3*ro2,r1*r4*ro3},{r1*r2*ro1,r2*r2,r2*r3*ro4,r2*r4*ro5},{r1*r3*ro2,r2*r3*ro4,r3*r3,r3*r4*ro6},{r1*r4*ro3,r2*r4*ro5,r3*r4*ro6,r4*r4} };
    arma::mat I = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    double help_s_1 = 1/(r1*r1*r2*r2*r3*r3*r4*r4*(1+2*ro4*ro5*ro6+2*ro1*ro2*ro4+2*ro1*ro3*ro5+2*ro2*ro3*ro6-ro1*ro1-ro2*ro2-ro3*ro3-ro4*ro4-ro5*ro5-ro6*ro6+ro1*ro1*ro6*ro6+ro2*ro2*ro5*ro5+ro3*ro3*ro4*ro4-2*ro1*ro3*ro4*ro6-2*ro1*ro2*ro5*ro6-2*ro2*ro3*ro4*ro5));
    arma::mat help_s_2 = {{r2*r2*r3*r3*r4*r4*(1+2*ro4*ro5*ro6-ro4*ro4-ro5*ro5-ro6*ro6),-r1*r2*r3*r3*r4*r4*(ro1+ro3*ro4*ro6+ro2*ro5*ro6-ro3*ro5-ro2*ro4-ro1*ro6*ro6),r1*r2*r2*r3*r4*r4*(ro1*ro4+ro3*ro6+ro2*ro5*ro5-ro3*ro4*ro5-ro2-ro1*ro5*ro6),-r1*r2*r2*r3*r3*r4*(ro1*ro4*ro6+ro3+ro2*ro4*ro5-ro3*ro4*ro4-ro2*ro6-ro1*ro5)},
    {-r1*r2*r3*r3*r4*r4*(ro1+ro3*ro4*ro6+ro2*ro5*ro6-ro3*ro5-ro2*ro4-ro1*ro6*ro6),r1*r1*r3*r3*r4*r4*(1+2*ro2*ro3*ro6-ro2*ro2-ro3*ro3-ro6*ro6),-r1*r1*r2*r3*r4*r4*(ro4+ro1*ro3*ro6+ro2*ro3*ro5-ro3*ro3*ro4-ro1*ro2-ro5*ro6),r1*r1*r2*r3*r3*r4*(ro4*ro6+ro1*ro3+ro2*ro2*ro5-ro2*ro3*ro4-ro1*ro2*ro6-ro5)},
    {r1*r2*r2*r3*r4*r4*(ro1*ro4+ro3*ro6+ro2*ro5*ro5-ro3*ro4*ro5-ro2-ro1*ro5*ro6),-r1*r1*r2*r3*r4*r4*(ro4+ro1*ro3*ro6+ro2*ro3*ro5-ro3*ro3*ro4-ro1*ro2-ro5*ro6),r1*r1*r2*r2*r4*r4*(1+2*ro1*ro3*ro5-ro1*ro1-ro3*ro3-ro5*ro5),-r1*r1*r2*r2*r3*r4*(ro6+ro1*ro3*ro4+ro1*ro2*ro5-ro2*ro3-ro1*ro1*ro6-ro4*ro5)},
    {-r1*r2*r2*r3*r3*r4*(ro1*ro4*ro6+ro3+ro2*ro4*ro5-ro3*ro4*ro4-ro2*ro6-ro1*ro5),r1*r1*r2*r3*r3*r4*(ro4*ro6+ro1*ro3+ro2*ro2*ro5-ro2*ro3*ro4-ro1*ro2*ro6-ro5),-r1*r1*r2*r2*r3*r4*(ro6+ro1*ro3*ro4+ro1*ro2*ro5-ro2*ro3-ro1*ro1*ro6-ro4*ro5),r1*r1*r2*r2*r3*r3*(1+2*ro1*ro2*ro4-ro1*ro1-ro2*ro2-ro4*ro4)}};
    arma::mat s_1 = help_s_2*help_s_1;
    for(int i = 0;i < m-1; i++){
      sigma_1.submat(span(4*i,4*i+3),span(4*i,4*i+3)) = s_1 + pow(v,(Times[i+1]-Times[i])) * s_1 * pow(v,(Times[i+1]-Times[i]));
      sigma_1.submat(span(4*i + 4,4*i + 7),span(4*i,4*i + 3)) = -s_1 * pow(v,(Times[i+1]-Times[i]));
      sigma_1.submat(span(4*i,4*i + 3),span(4*i + 4,4*i + 7)) = -pow(v,(Times[i+1]-Times[i])) * s_1;
    };
    sigma_1.submat(span(n-4,n-1),span(n-4,n-1)) = s_1;


    sigma_1_det = pow(help_s_1,m);

  }
  else{
    arma::vec phi = par(span(0,ranks-1));
    arma::vec r = par(span(ranks,2*ranks-1));
    arma::vec rho = par(span(2*ranks,((ranks+3)*ranks-2)/2));
    arma::mat v(ranks,ranks,fill::zeros);
    v.diag(0) = phi;
    arma::mat s(ranks,ranks);
    arma::mat x(ranks,ranks,fill::zeros);
    for(int i = 0;i<ranks;i++){
      s.col(i) = phi;
    }
    int help = 0;
    for(int i =0;i<ranks-1;i++){
        x(span(i+1,ranks-1),i) = rho(span(help,help + ranks -2-i));
        help = help + ranks - i - 1;
    }
    x = x.t() + x;
    arma::vec diag1(ranks,fill::ones);
    x.diag(0) = diag1;
    arma::mat sx = s*s.t()*x;
    arma::mat sx_1 = sx.i();
    double sx_1_det = det(sx_1);
    for(int i = 0;i<m-1;i++){
      sigma_1.submat(span(ranks*i,ranks*i+ranks-1),span(ranks*i,ranks*i+ranks-1)) = sx_1 + pow(v,(Times[i+1]-Times[i])) * sx_1 * pow(v,(Times[i+1]-Times[i]));
      sigma_1.submat(span(ranks*i + ranks,ranks*i + ranks+ranks-1),span(ranks*i,ranks*i + ranks-1)) = -sx_1 * pow(v,(Times[i+1]-Times[i]));
      sigma_1.submat(span(ranks*i,ranks*i + ranks-1),span(ranks*i + ranks,ranks*i + ranks+ranks-1)) = sigma_1.submat(span(ranks*i + ranks,ranks*i + ranks+ranks-1),span(ranks*i,ranks*i + ranks-1)).t();
    }
    sigma_1.submat(span(n-ranks,n-1),span(n-ranks,n-1)) = sx_1;
    sigma_1_det = pow(sx_1_det,m);
  }
  return Rcpp::List::create(Rcpp::Named("sigma_1") = wrap(sigma_1),
                              Rcpp::Named("det") = wrap(sigma_1_det));
}
//' H0 calculate mle
//' @param par numeric of mle par
//' @param N_phe input phenotype matrix
//' @param pheT numeric of times
//' @param ranks int, number of traits
//' @return  number of mle value
//'
//' @export
//[[Rcpp::export]]
double hfun_H0_mle (Rcpp::NumericVector par,arma::mat N_phe,Rcpp::NumericVector pheT,int ranks){
  Rcpp::NumericVector par_cov =wrap(as<arma::vec>(par)(span((3*ranks),((ranks*ranks+9*ranks-2)/2))));
  Rcpp::NumericVector par_cur =wrap(as<arma::vec>(par)(span(0,(3 * ranks - 1))));
  Rcpp::List cov_mat =hfun_get_invsigma(par_cov,pheT,ranks);
  arma::vec mu =hfun_get_mu(par_cur,pheT);

  int row = N_phe.n_rows;
  int col = N_phe.n_cols;
  arma::mat Mu(row,col,fill::zeros);
  for(int i = 0;i < row;i++){
    Mu.row(i) = mu.t();
  }
  arma::mat Y_delt = N_phe - Mu;
  double det = Rcpp::as<double>(cov_mat[1]);

  arma::mat Y_delt1 = cov_mat[0];
  arma::mat Y_delt2 = Y_delt * Y_delt1 % Y_delt;
  arma::vec Y_delt3(row);
  for(int i = 0;i<row;i++){
    Y_delt3(i) = sum(Y_delt2.row(i))/(-2);
  }
  Rcpp::NumericVector pv = Rcpp::as<Rcpp::NumericVector> (wrap(log(pow(1/(2*acos(double(-1))),col/2)*sqrt(det))+Y_delt3));
  double result = -sum(pv)/row;
  return result;
}
//' H1 calculate mle
//' @param par numeric of mle par
//' @param pheno list of pheno
//' @param Times numeric of times
//' @param ranks int, number of traits
//' @param ng int, number of geno type
//' @return  number of mle value
//'
//' @export
//[[Rcpp::export]]
double hfun_H1_mle (Rcpp::NumericVector par,Rcpp::List pheno,Rcpp::NumericVector Times,int ranks,int ng){
  Rcpp::NumericVector par_cov =wrap(as<arma::vec>(par)(span((3*ranks*ng),((ranks*ranks+3*ranks+6*ng*ranks-2)/2))));
  Rcpp::List cov_mat =hfun_get_invsigma(par_cov,Times,ranks);
  Rcpp::NumericVector Lh1(ng);
  int Row = 0;
  for(int i =0;i<ng;i++){
    Rcpp::NumericVector par_cur =wrap(as<arma::vec>(par)(span(i*3*ranks,((i+1)*3*ranks-1))));
    arma::vec mu = hfun_get_mu(par_cur,Times);
    arma::mat p = pheno[i];
    int row = p.n_rows;
    int col = p.n_cols;
    arma::mat Mu(row,col,fill::zeros);
    for(int i =0;i<row;i++){
      Mu.row(i) = mu.t();
    }
    arma::mat Y_delt = p - Mu;
    double det = Rcpp::as<double>(cov_mat[1]);
    arma::mat Y_delt1 = cov_mat[0];
    arma::mat Y_delt2 = Y_delt * Y_delt1 % Y_delt;
    arma::vec Y_delt3(row);
    for(int i =0;i<row;i++){
      Y_delt3(i) = sum(Y_delt2.row(i))/(-2);
    }
    Rcpp::NumericVector pv = Rcpp::as<Rcpp::NumericVector>(wrap(log(pow(1/(2*acos(double(-1))),col/2)*sqrt(det))+Y_delt3));
    Lh1[i] = -sum(pv);
    Row = Row + row;
  }
  return sum(Lh1)/Row;
}


/*** R




*/
