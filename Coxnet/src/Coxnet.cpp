// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;

//typedef Eigen::SparseVector<double> SpVec;
//typedef Eigen::SparseVector<double> SpVecf;
//typedef Eigen::SparseVector<int> SpVeci;
//typedef SpVec::InnerIterator InIterVec;
//typedef Eigen::SparseMatrix<double> SpMat;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  Eigen::VectorXd mX=X.colwise().mean(), sdX(X.cols()), sdXi;
  X.rowwise()-=mX.transpose();
  sdX=X.colwise().norm()/sqrt((double)X.rows());
  sdXi=1.0/sdX.array();
  X=X*sdXi.asDiagonal();
  return List::create(Named("x")=X, Named("sd")=sdX);
}
/*****  Soft-threshold  *****/
// [[Rcpp::export]]
double softC(double z,double lambda){
  if(z>lambda){
    return(z-lambda);
  }else if(z<-lambda){
    return(z+lambda);
  }else{
    return(0.0);
  }
}

/*****  Omega  *****/
// [[Rcpp::export]]
List OmegaC(Eigen::MatrixXd & Omega, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);
  
  //Omega.diagonal().setZero();
  Eigen::SparseMatrix<double> OmegaS=Omega.sparseView();
  
  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }
  
  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }
  
  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}

/*****  Sparse Omega  *****/
// [[Rcpp::export]]
List OmegaSC(Eigen::SparseMatrix<double> & OmegaS, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);
  
  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }
  
  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }
  
  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}





///////////////////////
/////  Non-local  /////
///////////////////////


/*****  Lambda path (max)  *****/
// [[Rcpp::export]]
double max_lambdaC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N,
Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, double alpha, Eigen::VectorXd wbeta, int N0){
  int i, j, q;
  double denS=N, c1=0.0;
  Eigen::VectorXd Li(N), lli(N);
  
  for(i=0;i<n;i++){
    c1+=(nevent1(i)/denS);denS-=nevent(i);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){lli(j)=tevent(j)-c1;}
  }
  Li=(lli.transpose()*X).cwiseAbs()/N0;
  Li=Li.array()/wbeta.array()/alpha;
  return Li.maxCoeff();
}

/*****  Derivatives of log-pl of eta (1st&2nd order),  ties  *****/
void dletaCm(Eigen::VectorXd& exb, Eigen::VectorXd& tevent, int& N, 
Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1,
int& n, Eigen::VectorXd& pl1, Eigen::VectorXd& pl2, int& ifast, int& itwo){
  int i, j, q, ipl2=0;
  double denSi, c1=0.0, c2=0.0;
  Eigen::VectorXd denS(n);
  
  if(ifast==0 || itwo==1)goto two;
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);c2+=(nevent1(i)/pow(denSi, 2));
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl1(j)=tevent(j)-exb(j)*c1;
      pl2(j)=exb(j)*(c1-exb(j)*c2);
      if(pl2(j)<=0.0)ipl2=1;
    }
  }
  if(ipl2==1){itwo=1;if(ifast==0){goto two;}}
  return;
  
  two:
  denSi=0.0;c1=0.0;c2=0.0;
  for(i=n-1;i>=0;--i){
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){denSi+=exb(j);}
    denS(i)=denSi;
  }
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denS(i));c2+=(nevent1(i)/pow(denS(i), 2));
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      pl1(j)=tevent(j)-exb(j)*c1;
      pl2(j)=exb(j)*(c1-exb(j)*c2);
    }
  }
}

/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent, 
Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo){
  int i, j, q, iSS=0;
  double ll=0.0, SSi;
  Eigen::VectorXd SS(n);

  if(ifast==0 || itwo==1)goto two;
  SSi=exb.sum();
  for(i=0;i<n;++i){
    if(SSi<=0.0)iSS=1;
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
    ll-=nevent1(i)*log(SSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
  }
  if(iSS==1){itwo=1;if(ifast==0){goto two;}}
  return(ll);
  
  two:
  ll=0.0;SSi=0.0;
  for(i=n-1;i>=0;--i){
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi+=exb(j);}
    SS(i)=SSi;
  }
  for(i=0;i<n;++i){
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){
      ll+=xb(j)-log(SS(i));
    }
  }
  return(ll);
}


/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvtrimC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, 
Eigen::MatrixXd XF, int NF, 
Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF, 
Eigen::MatrixXd X, int N, 
Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, 
int ifast, int itwo){
  int i, j;
  double lli, lfi;
  Eigen::VectorXd cv, xb=Eigen::VectorXd::Zero(N), xbF=Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd exb(N), exbF(NF);
  
  if(nn2>0){
    cv=Eigen::VectorXd::Zero(nn2);
    
    if(nn==0){
      exb=(xb.array()).exp();
      lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      exbF=(xbF.array()).exp();
      lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
      cv(0)=lfi-lli;
    }else{
      for(i=0;i<nn;i++){
        j=loco(i);
        xb+=X.col(j)*beta(i);exb=(xb.array()).exp();
        lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        xbF+=XF.col(j)*beta(i);exbF=(xbF.array()).exp();
        lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
        cv(i)=lfi-lli;
      }
    }
    if(nn2>nn){for(i=nn;i<nn2;i++){cv(i)=cv(nn-1);}}
  }else{
    cv=Eigen::VectorXd::Zero(1);
    
    exb=(xb.array()).exp();
    lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
    exbF=(xbF.array()).exp();
    lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
    cv(0)=lfi-lli;
  }
  
  return(cv);
}



/*****  Enet (L1+L2)  *****/
// [[Rcpp::export]]
List coxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, 
double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast){
  
  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQi=0.0;objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("nlambda")=il));
}

/*****  Enet (L1+L2)  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvcoxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, 
double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, 
int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){
  
  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQi=0.0;objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
    
    xbF=Eigen::VectorXd::Zero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*beta0(j);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}





/*****  Network (L1+La)  *****/
// [[Rcpp::export]]
List coxnetC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha, 
Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast){
  
  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        zi2=0.0;
        for(ij=0;ij<nadj(j);++ij){
          m=loc(ij, j)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
        }
        zi+=lambda2i*zi2;
        
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("nlambda")=il));
}

/*****  Network (L1+La)  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvcoxnetC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha, 
Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, 
int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){
  
  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        zi2=0.0;
        for(ij=0;ij<nadj(j);++ij){
          m=loc(ij, j)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
        }
        zi+=lambda2i*zi2;
        
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
    
    xbF=Eigen::VectorXd::Zero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*beta0(j);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}





///////////////////
/////  Local  /////
///////////////////

/*****  Max lambda path, local  *****/
// [[Rcpp::export]]
double max_loclambdaC(Eigen::MatrixXd X, Eigen::VectorXd tevent, Eigen::VectorXd Kh, Eigen::VectorXd Kh1,
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
int n, double alpha, Eigen::VectorXd wbeta, int N0){
  int i,j,q;
  double denS=Kh.sum(),c1=0.0;
  Eigen::VectorXd Li(N),lli(N);
  
  for(i=0;i<n;i++){
    c1+=(Kh1(i)/denS);
    for(j=loc1(i)-1,q=0;q<nevent(i);j++,q++){
      denS-=Kh(j);
      lli(j)=Kh(j)*tevent(j)-c1*Kh(j);
    }
  }
  Li=(lli.transpose()*X).cwiseAbs()/N0;
  Li=Li.array()/wbeta.array()/alpha;
  return Li.maxCoeff();
}

/*****  Derivatives of log-likelihood of eta (second order), local, ties  *****/
// [[Rcpp::export]]
Eigen::MatrixXd alocletaC(Eigen::ArrayXd eta,Eigen::VectorXd tevent,Eigen::VectorXd Kh,Eigen::VectorXd Kh1,int N,
Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,int n){
  int i,j,q;
  double denSi=0.0,c1=0.0,c2=0.0;
  Eigen::VectorXd eKh=eta.exp(),denS(n);
  Eigen::MatrixXd ll(N,2);
  
  eKh=eKh.cwiseProduct(Kh);
  for(i=n-1;i>=0;i--){
    for(j=loc1(i)-1,q=0;q<nevent(i);j++,q++){denSi+=eKh(j);}
    denS(i)=denSi;
  }
  
  for(i=0;i<n;i++){
    c1+=(Kh1(i)/denS(i));c2+=(Kh1(i)/pow(denS(i),2));
    for(j=loc1(i)-1,q=0;q<nevent(i);j++,q++){
      ll(j,1)=c2*pow(eKh(j),2)-c1*eKh(j);
      ll(j,0)=eta(j)-(Kh(j)*tevent(j)-c1*eKh(j))/ll(j,1);
    }
  }
  return(ll);
}

/*****  Log-likelihood of eta, local, ties  *****/
// [[Rcpp::export]]
double locletaC(Eigen::ArrayXd eta,Eigen::VectorXd Kh,
Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,int n){
  int i,j,q;
  double ll=0.0,SSi=0.0;
  Eigen::VectorXd eKh=eta.exp(),SS(n);
  
  eKh=eKh.cwiseProduct(Kh);
  for(i=n-1;i>=0;i--){
    for(j=loc1(i)-1,q=0;q<nevent(i);j++,q++){SSi+=eKh(j);}
    SS(i)=SSi;
  }
  for(i=0;i<n;i++){
    for(j=loc1(i)-1,q=0;q<nevent1(i);j++,q++){
      ll+=Kh(j)*(eta(j)-log(SS(i)));
    }
  }
  return(ll);
}


/*****  Objective function of Enet, local, ties  *****/
// [[Rcpp::export]]
double loceobjF(Eigen::VectorXd beta,Eigen::ArrayXd eta,Eigen::VectorXd lambda1,double lambda2,
Eigen::VectorXd Kh,Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,int n,int N0){
  double objQ;
  Eigen::SparseVector<double> sbeta=beta.sparseView();
  
  objQ=-locletaC(eta,Kh,nevent,nevent1,loc1,n)/N0;
  objQ+=sbeta.cwiseAbs().dot(lambda1);
  objQ+=sbeta.dot(sbeta)*lambda2/2.0;
  return(objQ);
}


/*****  Used for CV PL, local  *****/
// [[Rcpp::export]]
double loclbetaC(Eigen::VectorXd beta,Eigen::MatrixXd X,Eigen::VectorXd Kh,int N,
Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,int n){
  double ll;
  Eigen::VectorXd eta=Eigen::VectorXd::Zero(N);
  Eigen::SparseVector<double> sbeta=beta.sparseView();
  
  for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){eta+=X.col(i_.index())*i_.value();}
  ll=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
  return(ll);
}



/*****  Enet (L1+L2), local  *****/
// [[Rcpp::export]]
List locoxenetC(Eigen::MatrixXd X,Eigen::VectorXd tevent,
double alpha,Eigen::VectorXd lambda,int nlambda,Eigen::VectorXd wbeta,
Eigen::VectorXd Kh,Eigen::VectorXd Kh1,int N,Eigen::VectorXi nevent,
Eigen::VectorXi nevent1,Eigen::VectorXi loc1,
int n,int p,int N0,double thresh,double thresh2,int maxit){
  
  int i,it,il,io,ia=0,iai=0,II=0;
  double lambda2,zi,mi;
  int ii=0;
  Eigen::SparseVector<double> deno(p),sbeta;
  Eigen::SparseVector<int> ibeta(p);
  Eigen::SparseMatrix<double> Betas(p,nlambda);
  
  Eigen::MatrixXd ll(N,2);
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p),xb=Eigen::VectorXd::Zero(N),eta=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd lambdai(p),lambda1(p),Li(N),lli(N);
  Eigen::VectorXd locbeta=Eigen::VectorXd::Zero(nlambda),locbetai(maxit+1),obji(maxit+1);
  
  Eigen::VectorXi active(p),nzero=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd inactive=Eigen::VectorXd::Constant(p,1.0);
  Eigen::VectorXd::Index maxIndex;
  
  lambdai=lambda*alpha;
  
  for(il=0;il<nlambda;il++){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    
    io=0;
    sbeta=beta0.sparseView();
    locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
    obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+sbeta.dot(sbeta)*lambda2/2;
    
    it=0;
    while (1) {
      it++;
      
      ll=alocletaC(eta.array(),tevent,Kh,Kh1,N,nevent,nevent1,loc1,n);
      if((ll.col(1).array().abs()<thresh2).any()){II=4;break;}
      ll.col(1)=ll.col(1)/N0;
      lli=ll.col(1).cwiseProduct(ll.col(0));
      
      if(ia==0){
        Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X).cwiseAbs();
        Li=Li.array()/wbeta.array();
        do{
          mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
          if(mi>lambdai(il)){
            ibeta.insert(ia)=maxIndex;
            inactive(maxIndex)=0.0;ia++;
          }else{break;}
        }while(1);
      }      
      
      if(ia>0){
        do{
          ii++;
          active=ibeta.head(ia);
          for(i=0;i<ia;i++){
            deno.coeffRef(active(i))=lambda2-ll.col(1).dot(X.col(active(i)).cwiseAbs2());
            xb-=X.col(active(i))*beta0(active(i));
            zi=(xb.cwiseProduct(ll.col(1))-lli).dot(X.col(active(i)));
            beta0(active(i))=softC(zi,lambda1(active(i)))/deno.coeffRef(active(i));
            xb+=X.col(active(i))*beta0(active(i));
          }
          
          iai=ia;
          Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X).cwiseAbs();
          Li=Li.array()/wbeta.array();
          do{
            mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
            if(mi>lambdai(il)){
              ibeta.insert(ia)=maxIndex;
              inactive(maxIndex)=0.0;ia++;
            }else{break;}
          }while(1);
          if(iai==ia){break;} //ensure selection consistency
        }while(1);
      }
 
      eta=Eigen::VectorXd::Zero(N);
      for(i=0;i<ia;i++){eta+=X.col(active(i))*beta0(active(i));}
      
      io++;
      sbeta=beta0.sparseView();
      locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
      obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+sbeta.dot(sbeta)*lambda2/2;
      
      if (std::abs(locbetai(io-1)-locbetai(io))<std::abs(thresh*locbetai(io-1))){II=0;break;}
      if (std::abs(obji(io-1)-obji(io))<(thresh*obji(io-1))){II=0;break;}
      if (obji(io)!=obji(io)){II=2; goto exit;}
      if (it >= maxit) {II=1; goto exit;}
    } //while(it<maxit);
    
    //if(obji(io)!=obji(io))break;
    if((ll.col(1).array().abs()<thresh2).any()){II=3; goto exit;}
    
    locbeta(il)=locbetai(io);
    sbeta=beta0.sparseView();
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){Betas.insert(i_.index(),il)=i_.value();nzero(il)++;}
  }
  
  exit:
  Betas.makeCompressed();
  //if(II !=0 && il>0)--il;
  return(List::create(Named("Beta")=Betas,Named("ll")=locbeta,
  Named("flag")=II,Named("nzero")=nzero,Named("nlambda")=il));
}

/*****  Enet (L1+L2), local  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvlocoxenetC(Eigen::MatrixXd X,Eigen::VectorXd tevent,
double alpha,Eigen::VectorXd lambda,int nlambda,Eigen::VectorXd wbeta,
Eigen::VectorXd Kh,Eigen::VectorXd Kh1,int N,Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,
int n,int p,int N0,double thresh,double thresh2,int maxit,Eigen::MatrixXd XF,
Eigen::VectorXd KhF,int NF,Eigen::VectorXi neventF,Eigen::VectorXi nevent1F,Eigen::VectorXi loc1F,int nF){
  
  int i,it,il,io,ia=0,iai=0,II=0;
  double lambda2,zi,mi;
  int ii=0;
  Eigen::SparseVector<double> deno(p),sbeta;
  Eigen::SparseVector<int> ibeta(p);
  Eigen::SparseMatrix<double> Betas(p,nlambda);
  
  Eigen::MatrixXd ll(N,2);
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p),xb=Eigen::VectorXd::Zero(N),eta=Eigen::VectorXd::Zero(N),eta2(N);
  Eigen::VectorXd lambdai(p),lambda1(p),Li(N),lli(N);
  Eigen::VectorXd locbeta=Eigen::VectorXd::Zero(nlambda),locbetaF=Eigen::VectorXd::Zero(nlambda);
  Eigen::VectorXd locbetai(maxit+1),obji(maxit+1);
  
  Eigen::VectorXi active(p),nzero=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd inactive=Eigen::VectorXd::Constant(p,1.0);
  Eigen::VectorXd::Index maxIndex;
  
  lambdai=lambda*alpha;
  
  for(il=0;il<nlambda;il++){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    
    io=0;
    sbeta=beta0.sparseView();
    locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
    obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+sbeta.dot(sbeta)*lambda2/2;
    
    it=0;
    while (1) {
      it++;
      
      ll=alocletaC(eta.array(),tevent,Kh,Kh1,N,nevent,nevent1,loc1,n);
      if((ll.col(1).array().abs()<thresh2).any()){II=4;break;}
      ll.col(1)=ll.col(1)/N0;
      lli=ll.col(1).cwiseProduct(ll.col(0));
      
      if(ia==0){
        Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X).cwiseAbs();
        Li=Li.array()/wbeta.array();
        do{
          mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
          if(mi>lambdai(il)){
            ibeta.insert(ia)=maxIndex;
            inactive(maxIndex)=0.0;ia++;
          }else{break;}
        }while(1);
      }      
      
      if(ia>0){
        do{
          ii++;
          active=ibeta.head(ia);
          for(i=0;i<ia;i++){
            deno.coeffRef(active(i))=lambda2-ll.col(1).dot(X.col(active(i)).cwiseAbs2());
            xb-=X.col(active(i))*beta0(active(i));
            zi=(xb.cwiseProduct(ll.col(1))-lli).dot(X.col(active(i)));
            beta0(active(i))=softC(zi,lambda1(active(i)))/deno.coeffRef(active(i));
            xb+=X.col(active(i))*beta0(active(i));
          }
          
          iai=ia;
          Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X).cwiseAbs();
          Li=Li.array()/wbeta.array();
          do{
            mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
            if(mi>lambdai(il)){
              ibeta.insert(ia)=maxIndex;
              inactive(maxIndex)=0.0;ia++;
            }else{break;}
          }while(1);
          if(iai==ia){break;} //ensure selection consistency
        }while(1);
      }
 
      eta=Eigen::VectorXd::Zero(N);
      for(i=0;i<ia;i++){eta+=X.col(active(i))*beta0(active(i));}
      
      io++;
      sbeta=beta0.sparseView();
      locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
      obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+sbeta.dot(sbeta)*lambda2/2;
      
      if(std::abs(locbetai(io-1)-locbetai(io))<std::abs(thresh*locbetai(io-1))){II=0;break;}
      if(std::abs(obji(io-1)-obji(io))<(thresh*obji(io-1))){II=0;break;}
      if(obji(io)!=obji(io)){II=2;break;}
      if (it >= maxit) {II=1; goto exit;}
    } //while(it<maxit);
    
    //if(obji(io)!=obji(io))break;
    if((ll.col(1).array().abs()<thresh2).any()){II=3; goto exit;}
    
    locbeta(il)=locbetai(io);
    sbeta=beta0.sparseView();
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){Betas.insert(i_.index(),il)=i_.value();nzero(il)++;}
  
    eta2=Eigen::VectorXd::Zero(NF);
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){eta2+=XF.col(i_.index())*i_.value();}
    locbetaF(il)=locletaC(eta2.array(),KhF,neventF,nevent1F,loc1F,nF);
  }
  
  exit:
  Betas.makeCompressed();
  //if(II !=0 && il>0)--il;
  return(List::create(Named("Beta")=Betas,Named("ll")=locbeta,Named("lf")=locbetaF,
  Named("flag")=II,Named("nzero")=nzero,Named("nlambda")=il));
}





/*****  Network (L1+La), local  *****/
// [[Rcpp::export]]
List locoxnetC(Eigen::MatrixXd X,Eigen::VectorXd tevent,double alpha,Eigen::VectorXd lambda,int nlambda,
Eigen::VectorXd wbeta,Eigen::MatrixXd L,Eigen::MatrixXd Omega,
Eigen::VectorXd Kh,Eigen::VectorXd Kh1,int N,Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,
int n,int p,int N0,double thresh,double thresh2,int maxit){

  int i,it,il,io,ia=0,iai=0,II=0;
  double lambda2,zi,mi,objQi2;
  int ii=0;
  Eigen::SparseVector<double> deno(p),sbeta;
  Eigen::SparseVector<int> ibeta(p);
  Eigen::SparseMatrix<double> Betas(p,nlambda);
  
  Eigen::MatrixXd ll(N,2);
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p),xb=Eigen::VectorXd::Zero(N),eta=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd lambdai(p),lambda1(p),Li(N),lli(N);
  Eigen::VectorXd locbeta=Eigen::VectorXd::Zero(nlambda),locbetai(maxit+1),obji(maxit+1);
  
  Eigen::VectorXi active(p),nzero=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd inactive=Eigen::VectorXd::Constant(p,1.0);
  Eigen::VectorXd::Index maxIndex;
  
  lambdai=lambda*alpha;
  
  for(il=0;il<nlambda;il++){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    
    io=0;
    sbeta=beta0.sparseView();
    locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
    objQi2=beta0.transpose()*L*beta0;
    obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+objQi2*lambda2/2;
    
    it=0;
    while (1) {
      it++;
      
      ll=alocletaC(eta.array(),tevent,Kh,Kh1,N,nevent,nevent1,loc1,n);
      if((ll.col(1).array().abs()<thresh2).any()){II=4;break;}
      ll.col(1)=ll.col(1)/N0;
      lli=ll.col(1).cwiseProduct(ll.col(0));
      
      if(ia==0){
        Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X+lambda2*(beta0.transpose()*Omega)).cwiseAbs();
        Li=Li.array()/wbeta.array();
        do{
          mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
          if(mi>lambdai(il)){
            ibeta.insert(ia)=maxIndex;
            inactive(maxIndex)=0.0;ia++;
          }else{break;}
        }while(1);
      }
      
      if(ia>0){
        do{
          ii++;
          active=ibeta.head(ia);
          for(i=0;i<ia;i++){
            deno.coeffRef(active(i))=lambda2-ll.col(1).dot(X.col(active(i)).cwiseAbs2());
            xb-=X.col(active(i))*beta0(active(i));
            zi=(xb.cwiseProduct(ll.col(1))-lli).dot(X.col(active(i)))+lambda2*(beta0.dot(Omega.col(active(i))));
            beta0(active(i))=softC(zi,lambda1(active(i)))/deno.coeffRef(active(i));
            xb+=X.col(active(i))*beta0(active(i));
          }
          
          iai=ia;
          Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X+lambda2*(beta0.transpose()*Omega)).cwiseAbs();
          Li=Li.array()/wbeta.array();
          do{
            mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
            if(mi>lambdai(il)){
              ibeta.insert(ia)=maxIndex;
              inactive(maxIndex)=0.0;ia++;
            }else{break;}
          }while(1);
          if(iai==ia){break;} //ensure selection consistency
        }while(1);
      }
      
      eta=Eigen::VectorXd::Zero(N);
      for(i=0;i<ia;i++){eta+=X.col(active(i))*beta0(active(i));}
      
      io++;
      sbeta=beta0.sparseView();
      locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
      objQi2=beta0.transpose()*L*beta0;
      obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+objQi2*lambda2/2;
      
      if (std::abs(locbetai(io-1)-locbetai(io))<std::abs(thresh*locbetai(io-1))){II=0;break;}
      if (std::abs(obji(io-1)-obji(io))<(thresh*obji(io-1))){II=0;break;}
      if (obji(io)!=obji(io)){II=2;break;}
      if (it >= maxit) {II=1; goto exit;}
    }//while(it<maxit);
    
    //if(obji(io)!=obji(io))break;
    if((ll.col(1).array().abs()<thresh2).any()){II=3;break;}
    
    locbeta(il)=locbetai(io);
    sbeta=beta0.sparseView();
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){Betas.insert(i_.index(),il)=i_.value();nzero(il)++;}
  }
  exit:
  Betas.makeCompressed();
  //if(II !=0 && il>0)--il;
  return(List::create(Named("Beta")=Betas,Named("ll")=locbeta,
  Named("flag")=II,Named("nzero")=nzero,Named("nlambda")=il));
}

/*****  Network (L1+La), local  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvlocoxnetC(Eigen::MatrixXd X,Eigen::VectorXd tevent,double alpha,Eigen::VectorXd lambda,int nlambda,
Eigen::VectorXd wbeta,Eigen::MatrixXd L,Eigen::MatrixXd Omega,
Eigen::VectorXd Kh,Eigen::VectorXd Kh1,int N,Eigen::VectorXi nevent,Eigen::VectorXi nevent1,Eigen::VectorXi loc1,
int n,int p,int N0,double thresh,double thresh2,int maxit,Eigen::MatrixXd XF,
Eigen::VectorXd KhF,int NF,Eigen::VectorXi neventF,Eigen::VectorXi nevent1F,Eigen::VectorXi loc1F,int nF){
  
  int i,it,il,io,ia=0,iai=0,II=0;
  double lambda2,zi,mi,objQi2;
  int ii=0;
  Eigen::SparseVector<double> deno(p),sbeta;
  Eigen::SparseVector<int> ibeta(p);
  Eigen::SparseMatrix<double> Betas(p,nlambda);
  
  Eigen::MatrixXd ll(N,2);
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p),xb=Eigen::VectorXd::Zero(N),eta=Eigen::VectorXd::Zero(N),eta2(N);
  Eigen::VectorXd lambdai(p),lambda1(p),Li(N),lli(N);
  Eigen::VectorXd locbeta=Eigen::VectorXd::Zero(nlambda),locbetaF=Eigen::VectorXd::Zero(nlambda);
  Eigen::VectorXd locbetai(maxit+1),obji(maxit+1);
  
  Eigen::VectorXi active(p),nzero=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd inactive=Eigen::VectorXd::Constant(p,1.0);
  Eigen::VectorXd::Index maxIndex;
  
  lambdai=lambda*alpha;
  
  for(il=0;il<nlambda;il++){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    
    io=0;
    sbeta=beta0.sparseView();
    locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
    objQi2=beta0.transpose()*L*beta0;
    obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+objQi2*lambda2/2;
    
    it=0;
    while (1) {
      it++;
      
      ll=alocletaC(eta.array(),tevent,Kh,Kh1,N,nevent,nevent1,loc1,n);
      if((ll.col(1).array().abs()<thresh2).any()){II=4;break;}
      ll.col(1)=ll.col(1)/N0;
      lli=ll.col(1).cwiseProduct(ll.col(0));
      
      if(ia==0){
        Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X+lambda2*(beta0.transpose()*Omega)).cwiseAbs();
        Li=Li.array()/wbeta.array();
        do{
          mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
          if(mi>lambdai(il)){
            ibeta.insert(ia)=maxIndex;
            inactive(maxIndex)=0.0;ia++;
          }else{break;}
        }while(1);
      }
      
      if(ia>0){
        do{
          ii++;
          active=ibeta.head(ia);
          for(i=0;i<ia;i++){
            deno.coeffRef(active(i))=lambda2-ll.col(1).dot(X.col(active(i)).cwiseAbs2());
            xb-=X.col(active(i))*beta0(active(i));
            zi=(xb.cwiseProduct(ll.col(1))-lli).dot(X.col(active(i)))+lambda2*(beta0.dot(Omega.col(active(i))));
            beta0(active(i))=softC(zi,lambda1(active(i)))/deno.coeffRef(active(i));
            xb+=X.col(active(i))*beta0(active(i));
          }
          
          iai=ia;
          Li=((xb.cwiseProduct(ll.col(1))-lli).transpose()*X+lambda2*(beta0.transpose()*Omega)).cwiseAbs();
          Li=Li.array()/wbeta.array();
          do{
            mi=Li.cwiseProduct(inactive).maxCoeff(&maxIndex);
            if(mi>lambdai(il)){
              ibeta.insert(ia)=maxIndex;
              inactive(maxIndex)=0.0;ia++;
            }else{break;}
          }while(1);
          if(iai==ia){break;} //ensure selection consistency
        }while(1);
      }
      
      eta=Eigen::VectorXd::Zero(N);
      for(i=0;i<ia;i++){eta+=X.col(active(i))*beta0(active(i));}
      
      io++;
      sbeta=beta0.sparseView();
      locbetai(io)=locletaC(eta.array(),Kh,nevent,nevent1,loc1,n);
      objQi2=beta0.transpose()*L*beta0;
      obji(io)=-locbetai(io)/N0+sbeta.cwiseAbs().dot(lambda1)+objQi2*lambda2/2;
      
      if (std::abs(locbetai(io-1)-locbetai(io))<std::abs(thresh*locbetai(io-1))){II=0;break;}
      if (std::abs(obji(io-1)-obji(io))<(thresh*obji(io-1))){II=0;break;}
      if (obji(io)!=obji(io)){II=2; goto exit;}
      if (it >= maxit) {II=1; goto exit;}
    }//while(it<maxit);
    
    //if(obji(io)!=obji(io))break;
    if((ll.col(1).array().abs()<thresh2).any()){II=3;break;}
    
    locbeta(il)=locbetai(io);
    sbeta=beta0.sparseView();
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){Betas.insert(i_.index(),il)=i_.value();nzero(il)++;}
    
    eta2=Eigen::VectorXd::Zero(NF);
    for(Eigen::SparseVector<double>::InnerIterator i_(sbeta);i_;++i_){eta2+=XF.col(i_.index())*i_.value();}
    locbetaF(il)=locletaC(eta2.array(),KhF,neventF,nevent1F,loc1F,nF);
  }
  
  exit:
  Betas.makeCompressed();
  //if(II !=0 && il>0)--il;
  return(List::create(Named("Beta")=Betas,Named("ll")=locbeta,Named("lf")=locbetaF,
  Named("flag")=II,Named("nzero")=nzero,Named("nlambda")=il));
}




