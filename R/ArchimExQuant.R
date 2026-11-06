ArchimExQuant=function(tau,X,Y,x,k="kopt",copula,ci.level=0.95,method="Gaussian",smoothing=FALSE,nboot=1000){

  n=length(Y)

  X=as.matrix(X)

  if(dim(X)[1]!=length(Y)){
    stop("X and Y must have the same length !")
  }

  if(dim(X)[2]!=length(x)){
    stop("X and x must have the same length !")
  }

  if(!is.logical(smoothing)){
    stop("smoothing must be boolean.")
  }

  if(!is.numeric(nboot) || nboot<=0){
    stop("nboot must be a strictly positive number")
  }

  if (k!="kopt" && k > n - 1 || k < 1) {
    stop("k must be between 1 and n-1.")
  }

  if (ci.level >= 1 || ci.level <= 0) {
    stop("ci.level must be (strictly) between 0 and 1.")
  }

  if (tau >= 1 || tau <= 0) {
    stop("tau must be (strictly) between 0 and 1.")
  }

  if (copula != "Clayton" && copula != "Frank" && copula != "Gumbel" && copula != "Joe") {
    stop("copula may be either Clayton, Frank, Gumbel or Joe.")
  }

  if (method != "Gaussian" && method != "Gaussian_Unified" && method != "Bootstrap") {
    stop("method may be either Gaussian, Gaussian_Unified or Bootstrap.")
  }

  X=as.matrix(X)

  d=dim(X)[2]+1

  mop_br=mop(Y, 1:(n - 1), 0, method = "RBMOP")
  rhohat=mop_br$rho
  bhat=mop_br$beta

  if(k[1]=="kopt"){
    k=trunc(((1-rhohat)^2/(-2*rhohat*bhat^2))^(1/(1-2*rhohat))*n^(-2*rhohat/(1-2*rhohat)))
  }

  if(smoothing==FALSE){
    tauhat=apply(X<=matrix(rep(x,n),n,d-1,byrow=T),2,mean)
  }

  if(smoothing==TRUE){
    hn=apply(as.matrix(X),2,sd)*(4/3)^(1/5)*n^(-1/5)
    tauhat=apply(pnorm((matrix(rep(x,n),n,d-1,byrow=T)-X)/matrix(rep(hn,n),n,d-1,byrow=T)),2,mean)
  }

  kendall_tau=(sum(cor(cbind(Y,X),method="kendall"))-d)/(d*(d-1))

  if(copula=="Clayton"){

    deltahat=2*kendall_tau/(1-kendall_tau)

    phi=function(t){
      return((t^(-deltahat)-1)/deltahat)
    }

    phiinv=function(t){
      return((1+deltahat*t)^(-1/deltahat))
    }

    psi=function(t){
      Cd=(-1)^(d-1)*prod(1+(1:(d-1)-1)*deltahat)
      return(Cd*(1+deltahat*t)^(-1/deltahat-(d-1)))
    }

    psid=function(t){
      Cd=(-1)^(d)*prod(1+(1:(d)-1)*deltahat)
      return(Cd*(1+deltahat*t)^(-1/deltahat-(d)))
    }

    psiinv=function(t){
      Cd=(-1)^(d-1)*prod(1+(1:(d-1)-1)*deltahat)
      return(((t/Cd)^(-deltahat/(1+(d-1)*deltahat))-1)/deltahat)
    }

  }

  if(copula=="Frank"){

    f=function(x){
      return(1-4/x*(1-debye1(x))-kendall_tau)
    }

    deltahat=uniroot(f,c(1e-10,100))$root

    phi=function(t){
      return(-log((exp(-deltahat*t)-1)/(exp(-deltahat)-1)))
    }

    phiinv=function(t){
      return(-1/deltahat*log(1+(exp(-deltahat)-1)*exp(-t)))
    }

    psi=function(t) -1/deltahat*log(1+(exp(-deltahat)-1)*exp(-t))

    id=2

    while(id<=d){

      psi=Deriv(psi)

      id=id+1
    }

    psid=Deriv(psi)
    psid1=Deriv(psid)

    psiinv=function(t){
      xk=(tau-1)*psi(omegahat)/psid(omegahat)-(1-tau)^2*psi(omegahat)^2/2*psid1(omegahat)/psid(omegahat)^3
      i=1
      while(i<=100){
        xk=xk-(psi(xk)-t)/psid(xk)
        i=i+1
      }
      return(xk)
    }

  }

  if(copula=="Gumbel"){

    deltahat=1/(1-kendall_tau)

    phi=function(t){
      return((-log(t))^deltahat)
    }

    phiinv=function(t){
      return(exp(-t^(1/deltahat)))
    }

    psi=function(t) exp(-t^(1/deltahat))

    id=2

    while(id<=d){

      psi=Deriv(psi)

      id=id+1
    }

    psid=Deriv(psi)
    psid1=Deriv(psid)

    psiinv=function(t){
      xk=(tau-1)*psi(omegahat)/psid(omegahat)-(1-tau)^2*psi(omegahat)^2/2*psid1(omegahat)/psid(omegahat)^3
      i=1
      while(i<=100){
        xk=xk-(psi(xk)-t)/psid(xk)
        i=i+1
      }
      return(xk)
    }

  }

  if(copula=="Joe"){

    f=function(x){
      return(1+2/(2-x)*(digamma(2)-digamma(1+2/x))-kendall_tau)
    }

    deltahat=uniroot(f,c(1,1000))$root

    phi=function(t){
      return(-log(1-(1-t)^deltahat))
    }

    phiinv=function(t){
      return(1-(1-exp(-t))^(1/deltahat))
    }

    psi=function(t) 1-(1-exp(-t))^(1/deltahat)

    id=2

    while(id<=d){

      psi=Deriv(psi)

      id=id+1
    }

    psid=Deriv(psi)
    psid1=Deriv(psid)

    psiinv=function(t){
      xk=(tau-1)*psi(omegahat)/psid(omegahat)-(1-tau)^2*psi(omegahat)^2/2*psid1(omegahat)/psid(omegahat)^3
      i=1
      while(i<=100){
        xk=xk-(psi(xk)-t)/psid(xk)
        i=i+1
      }
      return(xk)
    }

  }

  omegahat=sum(phi(tauhat))
  tauprime=phiinv(psiinv(tau*psi(omegahat))-omegahat)

  if(copula=="Clayton" || copula=="Frank"){
    M=(-psi(omegahat)/psid(omegahat))^(-1)
  }

  if(copula=="Gumbel" || copula=="Joe"){
    M=(-psi(omegahat)/psid(omegahat))^(-1/deltahat)
  }

  if(method=="Gaussian"){
    k=max(k,trunc(n*(1-phiinv(1-tau)))+1)
    l=k/(n*(1-phiinv(1-tau)))
  }

  if(method=="Gaussian_Unified"){
    k=max(k,trunc(n*(1-phiinv(-(1-tau)*psi(omegahat)/psid(omegahat))))+1)
    lbar=k/(n*(1-phiinv(-(1-tau)*psi(omegahat)/psid(omegahat))))
  }

  if(method=="Bootstrap"){
    k=max(k,trunc(n*(1-tauprime))+1)
  }

  gammahat = mop(Y, k, 0, method = "MOP")$EVI
  quanthat=as.vector(quantile(Y,1-k/n))*(k/(n*(1-tauprime)))^as.vector(gammahat)

  if(method=="Gaussian"){
    quantup=as.vector(quanthat*exp(log(l)/sqrt(k)*sqrt((1+log(M)/log(l))^2+1/log(l)^2)*gammahat*qnorm((1+ci.level)/2)))
    quantdown=as.vector(quanthat*exp(-log(l)/sqrt(k)*sqrt((1+log(M)/log(l))^2+1/log(l)^2)*gammahat*qnorm((1+ci.level)/2)))
  }

  if(method=="Gaussian_Unified"){
    quantup=as.vector(quanthat*exp(log(lbar)/sqrt(k)*sqrt(1+1/log(lbar)^2)*gammahat*qnorm((1+ci.level)/2)))
    quantdown=as.vector(quanthat*exp(-log(lbar)/sqrt(k)*sqrt(1+1/log(lbar)^2)*gammahat*qnorm((1+ci.level)/2)))
  }

  if(method=="Bootstrap"){

    set.seed(2025)
    quantboot=rep(0,nboot)

    if(copula=="Clayton"){

      phi_boot=function(t){
        return((t^(-deltahat_boot)-1)/deltahat_boot)
      }

      phiinv_boot=function(t){
        return((1+deltahat_boot*t)^(-1/deltahat_boot))
      }

      psi_boot=function(t){
        Cd=(-1)^(d-1)*prod(1+(1:(d-1)-1)*deltahat_boot)
        return(Cd*(1+deltahat_boot*t)^(-1/deltahat_boot-(d-1)))
      }

      psiinv_boot=function(t){
        Cd=(-1)^(d-1)*prod(1+(1:(d-1)-1)*deltahat_boot)
        return(((t/Cd)^(-deltahat_boot/(1+(d-1)*deltahat_boot))-1)/deltahat_boot)
      }

      j=1

      while(j<=nboot){

        Z_boot=rCopula(n,claytonCopula(param=deltahat,dim=d))
        Y_boot=quantile(Y,Z_boot[,1])*(Z_boot[,1]<=1-k/n)+quantile(Y,1-k/n,type=1)*(n*(1-Z_boot[,1])/k)^(-as.vector(gammahat))*(Z_boot[,1]>1-k/n)

        X_boot=matrix(0,n,d-1)

        i=2
        while(i<=d){

          X_boot[,i-1]=quantile(X[,i-1],Z_boot[,i])

          i=i+1
        }

        kendall_tau_boot=(sum(cor(Z_boot,method="kendall"))-d)/(d*(d-1))

        deltahat_boot=2*kendall_tau_boot/(1-kendall_tau_boot)

        if(smoothing==FALSE){
          tauhat_boot=apply(X_boot<=matrix(rep(x,n),n,length(x),byrow=T),2,mean)
        }

        if(smoothing==TRUE){
          hn_boot=apply(as.matrix(X_boot),2,sd)*(4/3)^(1/5)*n^(-1/5)
          tauhat_boot=apply(pnorm((matrix(rep(x,n),n,d-1,byrow=T)-X_boot)/matrix(rep(hn_boot,n),n,d-1,byrow=T)),2,mean)
        }

        omegahat_boot=sum(phi_boot(tauhat_boot))
        mopest_boot=mop(Y_boot,k=k,p=0,method="MOP")
        gammahat_boot=mopest_boot$EVI

        tauprime_boot=phiinv_boot(psiinv_boot(tau*psi_boot(omegahat_boot))-omegahat_boot)

        mop_boot=mop(Y_boot, 1:(n - 1), 0, method = "RBMOP")
        rho_boot=mop_boot$rho
        b_boot=mop_boot$beta

        k_boot=min(trunc(((1-rho_boot)^2/(-2*rho_boot*b_boot^2))^(1/(1-2*rho_boot))*n^(-2*rho_boot/(1-2*rho_boot))),n-1)
        k_boot=max(k_boot,trunc(n*(1-tauprime_boot))+1)

        quantboot[j]=quantile(Y_boot,1-k_boot/n)*(k_boot/(n*(1-tauprime_boot)))^gammahat_boot

        j=j+1
      }

    }
    if(copula=="Frank"){

      phi_boot=function(t){
        return(-log((exp(-deltahat_boot*t)-1)/(exp(-deltahat_boot)-1)))
      }

      phiinv_boot=function(t){
        return(-1/deltahat_boot*log(1+(exp(-deltahat_boot)-1)*exp(-t)))
      }

      psi_boot=function(t) -1/deltahat_boot*log(1+(exp(-deltahat_boot)-1)*exp(-t))

      id=2

      while(id<=d){

        psi_boot=Deriv(psi_boot)

        id=id+1
      }

      psid_boot=Deriv(psi_boot)
      psid1_boot=Deriv(psid_boot)

      psiinv_boot=function(t){
        xk=(tau-1)*psi_boot(omegahat_boot)/psid_boot(omegahat_boot)-(1-tau)^2*psi_boot(omegahat_boot)^2/2*psid1_boot(omegahat_boot)/psid_boot(omegahat_boot)^3
        i=1
        while(i<=100){
          xk=xk-(psi_boot(xk)-t)/psid_boot(xk)
          i=i+1
        }
        return(xk)
      }

      j=1

      while(j<=nboot){

        Z_boot=rCopula(n,frankCopula(param=deltahat,dim=d))
        Y_boot=quantile(Y,Z_boot[,1])*(Z_boot[,1]<=1-k/n)+quantile(Y,1-k/n,type=1)*(n*(1-Z_boot[,1])/k)^(-as.vector(gammahat))*(Z_boot[,1]>1-k/n)

        X_boot=matrix(0,n,d-1)

        i=2
        while(i<=d){

          X_boot[,i-1]=quantile(X[,i-1],Z_boot[,i])

          i=i+1
        }

        kendall_tau_boot=(sum(cor(Z_boot,method="kendall"))-d)/(d*(d-1))

        f_boot=function(x){
          return(1-4/x*(1-debye1(x))-kendall_tau_boot)
        }

        deltahat_boot=uniroot(f_boot,c(1e-10,100))$root

        if(smoothing==FALSE){
          tauhat_boot=apply(X_boot<=matrix(rep(x,n),n,length(x),byrow=T),2,mean)
        }

        if(smoothing==TRUE){
          hn_boot=apply(as.matrix(X_boot),2,sd)*(4/3)^(1/5)*n^(-1/5)
          tauhat_boot=apply(pnorm((matrix(rep(x,n),n,d-1,byrow=T)-X_boot)/matrix(rep(hn_boot,n),n,d-1,byrow=T)),2,mean)
        }

        omegahat_boot=sum(phi_boot(tauhat_boot))
        mopest_boot=mop(Y_boot,k=k,p=0,method="MOP")
        gammahat_boot=mopest_boot$EVI

        tauprime_boot=phiinv_boot(psiinv_boot(tau*psi_boot(omegahat_boot))-omegahat_boot)

        mop_boot=mop(Y_boot, 1:(n - 1), 0, method = "RBMOP")
        rho_boot=mop_boot$rho
        b_boot=mop_boot$beta

        k_boot=min(trunc(((1-rho_boot)^2/(-2*rho_boot*b_boot^2))^(1/(1-2*rho_boot))*n^(-2*rho_boot/(1-2*rho_boot))),n-1)
        k_boot=max(k_boot,trunc(n*(1-tauprime_boot))+1)

        quantboot[j]=quantile(Y_boot,1-k_boot/n)*(k_boot/(n*(1-tauprime_boot)))^gammahat_boot

        j=j+1
      }

    }
    if(copula=="Gumbel"){

      phi_boot=function(t){
        return((-log(t))^deltahat_boot)
      }

      phiinv_boot=function(t){
        return(exp(-t^(1/deltahat_boot)))
      }

      psi_boot=function(t) exp(-t^(1/deltahat_boot))

      id=2

      while(id<=d){

        psi_boot=Deriv(psi_boot)

        id=id+1
      }

      psid_boot=Deriv(psi_boot)
      psid1_boot=Deriv(psid_boot)

      psiinv_boot=function(t){
        xk=(tau-1)*psi_boot(omegahat_boot)/psid_boot(omegahat_boot)-(1-tau)^2*psi_boot(omegahat_boot)^2/2*psid1_boot(omegahat_boot)/psid_boot(omegahat_boot)^3
        i=1
        while(i<=100){
          xk=xk-(psi_boot(xk)-t)/psid_boot(xk)
          i=i+1
        }
        return(xk)
      }

      j=1

      while(j<=nboot){

        Z_boot=rCopula(n,gumbelCopula(param=deltahat,dim=d))
        Y_boot=quantile(Y,Z_boot[,1])*(Z_boot[,1]<=1-k/n)+quantile(Y,1-k/n,type=1)*(n*(1-Z_boot[,1])/k)^(-as.vector(gammahat))*(Z_boot[,1]>1-k/n)

        X_boot=matrix(0,n,d-1)

        i=2
        while(i<=d){

          X_boot[,i-1]=quantile(X[,i-1],Z_boot[,i])

          i=i+1
        }

        kendall_tau_boot=(sum(cor(Z_boot,method="kendall"))-d)/(d*(d-1))

        deltahat_boot=1/(1-kendall_tau_boot)

        if(smoothing==FALSE){
          tauhat_boot=apply(X_boot<=matrix(rep(x,n),n,length(x),byrow=T),2,mean)
        }

        if(smoothing==TRUE){
          hn_boot=apply(as.matrix(X_boot),2,sd)*(4/3)^(1/5)*n^(-1/5)
          tauhat_boot=apply(pnorm((matrix(rep(x,n),n,d-1,byrow=T)-X_boot)/matrix(rep(hn_boot,n),n,d-1,byrow=T)),2,mean)
        }

        omegahat_boot=sum(phi_boot(tauhat_boot))
        mopest_boot=mop(Y_boot,k=k,p=0,method="MOP")
        gammahat_boot=mopest_boot$EVI

        tauprime_boot=phiinv_boot(psiinv_boot(tau*psi_boot(omegahat_boot))-omegahat_boot)

        mop_boot=mop(Y_boot, 1:(n - 1), 0, method = "RBMOP")
        rho_boot=mop_boot$rho
        b_boot=mop_boot$beta

        k_boot=min(trunc(((1-rho_boot)^2/(-2*rho_boot*b_boot^2))^(1/(1-2*rho_boot))*n^(-2*rho_boot/(1-2*rho_boot))),n-1)
        k_boot=max(k_boot,trunc(n*(1-tauprime_boot))+1)

        quantboot[j]=quantile(Y_boot,1-k_boot/n)*(k_boot/(n*(1-tauprime_boot)))^gammahat_boot

        j=j+1
      }

    }
    if(copula=="Joe"){

      phi_boot=function(t){
        return(-log(1-(1-t)^deltahat_boot))
      }

      phiinv_boot=function(t){
        return(1-(1-exp(-t))^(1/deltahat_boot))
      }

      psi_boot=function(t) 1-(1-exp(-t))^(1/deltahat_boot)

      id=2

      while(id<=d){

        psi_boot=Deriv(psi_boot)

        id=id+1
      }

      psid_boot=Deriv(psi_boot)
      psid1_boot=Deriv(psid_boot)

      psiinv_boot=function(t){
        xk=(tau-1)*psi_boot(omegahat_boot)/psid_boot(omegahat_boot)-(1-tau)^2*psi_boot(omegahat_boot)^2/2*psid1_boot(omegahat_boot)/psid_boot(omegahat_boot)^3
        i=1
        while(i<=100){
          xk=xk-(psi_boot(xk)-t)/psid_boot(xk)
          i=i+1
        }
        return(xk)
      }

      j=1

      while(j<=nboot){

        Z_boot=rCopula(n,joeCopula(param=deltahat,dim=d))
        Y_boot=quantile(Y,Z_boot[,1])*(Z_boot[,1]<=1-k/n)+quantile(Y,1-k/n,type=1)*(n*(1-Z_boot[,1])/k)^(-as.vector(gammahat))*(Z_boot[,1]>1-k/n)

        X_boot=matrix(0,n,d-1)

        i=2
        while(i<=d){

          X_boot[,i-1]=quantile(X[,i-1],Z_boot[,i])

          i=i+1
        }

        kendall_tau_boot=(sum(cor(Z_boot,method="kendall"))-d)/(d*(d-1))

        f_boot=function(x){
          return(1+2/(2-x)*(digamma(2)-digamma(1+2/x))-kendall_tau_boot)
        }

        deltahat_boot=uniroot(f_boot,c(1,1000))$root

        if(smoothing==FALSE){
          tauhat_boot=apply(X_boot<=matrix(rep(x,n),n,length(x),byrow=T),2,mean)
        }

        if(smoothing==TRUE){
          hn_boot=apply(as.matrix(X_boot),2,sd)*(4/3)^(1/5)*n^(-1/5)
          tauhat_boot=apply(pnorm((matrix(rep(x,n),n,d-1,byrow=T)-X_boot)/matrix(rep(hn_boot,n),n,d-1,byrow=T)),2,mean)
        }

        omegahat_boot=sum(phi_boot(tauhat_boot))
        mopest_boot=mop(Y_boot,k=k,p=0,method="MOP")
        gammahat_boot=mopest_boot$EVI

        tauprime_boot=phiinv_boot(psiinv_boot(tau*psi_boot(omegahat_boot))-omegahat_boot)

        mop_boot=mop(Y_boot, 1:(n - 1), 0, method = "RBMOP")
        rho_boot=mop_boot$rho
        b_boot=mop_boot$beta

        k_boot=min(trunc(((1-rho_boot)^2/(-2*rho_boot*b_boot^2))^(1/(1-2*rho_boot))*n^(-2*rho_boot/(1-2*rho_boot))),n-1)
        k_boot=max(k_boot,trunc(n*(1-tauprime_boot))+1)

        quantboot[j]=quantile(Y_boot,1-k_boot/n)*(k_boot/(n*(1-tauprime_boot)))^gammahat_boot

        j=j+1
      }

    }

    quantup=as.vector(quantile(quantboot,1-(1-ci.level)/2,na.rm=T))
    quantdown=as.vector(quantile(quantboot,(1-ci.level)/2,na.rm=T))

  }

  return(list(Parameter_estimate=deltahat,k=k,Lower_bound=quantdown,Point_estimate=quanthat,Upper_bound=quantup))

}
