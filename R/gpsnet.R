get_lap<-function(A,cor_m){
  
  cc<-diag(colSums(A))
  n<-length(cc[,1])
  L1<-diag(1,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n) {
      if(A[i,j]!=0) 
      {L1[i,j]=-1/sqrt(cc[i,i]*cc[j,j])}
    }
  }
  L1=L1*abs(cor_m)
  
  Le<-forceSymmetric(L1)
  Le
}

sim_1<-function(N,gp,gn,vv)
{
  p<-gp*gn
  Lc<-list()
  Al<-list()
  le<-list()
  datax<-c()
  K<-c()
  
  for (i in 1:gn) {
    l1<-huge.generator(n=N,d=gp,prob = vv,vis=TRUE)
    x1<-l1$data
    Lc[[i]]<-l1$sigma
    Al[[i]]<-as.matrix(l1$theta)
    l2<-huge(x1)
    lee<-l2[[5]]
    datax<-cbind(datax,x1)
    cr<-cor(x1)
    Ae<-l2$path[[5]]
    Ae<-as.matrix(Ae)
    le[[i]]<-get_lap(Ae,cr)
    linv<-solve(le[[i]])
    uu<-x1%*%linv%*%t(x1)
    K[[i]]<-uu
    
  }
  list_x<-list(sim_x=datax, la_e=le,t_cv=Lc,tr_p=Al,Kernel=K )
  list_x
  
}
sim_sv_g<-function(datax,c_rate=0.3,q,gp,gn,bs,tg)
{
  N=length(datax[,1])
  p=length(datax[1,])
  x=as.matrix(datax)
  vu=1
  lambda=1
  beta_g<-rep(0,gp)
  beta_g[1:q]<-bs
  beta0<-rep(beta_g,tg)
  beta<-rep(0,gp*gn)
  beta[1:length(beta0)]<-beta0
  mu<-exp(drop(x%*%beta))*2
  real_time<--(log(runif(N)))/((mu)^(1/vu))
  if(c_rate<1)
  {lam_c<-(lambda*c_rate/(1-c_rate))^vu}
  
  mu_c<-exp(drop(x %*% beta))*lam_c
  fl_time<--(log(runif(N)))/((mu_c)^(1/vu))
  status<-as.numeric(real_time<=fl_time)
  if(c_rate>=1) status=rep(0,size)
  time<-pmin(real_time,fl_time)
  sim_dat<-data.frame(time,status)
  x<-as.data.frame(datax)
  re<-list(datax=x,datay=sim_dat)
  re
}

sim_v1<-function(N,gp,gn,q,bs,tg,vv,c_rate)
{
  lsim=sim_1(N,gp,gn,vv)
  dd<-sim_sv_g(datax=lsim$sim_x,q = q,gp=gp,gn=gn,bs=bs,tg=tg)
  list_sim<-list(datax=dd$datax,datay=dd$datay,net_e=lsim$la_e,true_path=lsim$tr_p, true_cov=lsim$t_cv,Ker=lsim$Kernel)  
  
}
###data split
da_gp<-function(datax,gn,gp)
{
  da<-list()
  for (i in 1:gn) {
    idx1<-1+gp*(i-1)
    idx2<-i*gp
    da[[i]]<-datax[,idx1:idx2]
  }
  da
}

####################################kerenl
find_sz<-function(a,b,c)
{
  
  id_b<-which(b>0)
  ua<-a[id_b]
  ub<-b[id_b]
  uc<-c[id_b]
  uab<-ua+ub
  uc_ab<-uc-uab
  if(min(uc_ab)>0) {sz=1}
  else
  {u<-(uc-ua)/ub
  
  um<-min(u)
  sz<-0.999*um}
  
  
  sz
  
}


fit_mk_cox<-function(K,datax,datay,gp,gn,C,lambda)
{
  nd<-cox_arr(datax,datay)
  N=length(datax[,1])
  dg<-da_gp(nd[[1]],gn,gp)
  id_new<-nd[[3]]
  #kernel adjust with re-ordered data
  for (m in 1:gn) {
    kk=K[[i]]
    kk=kk[id_new,]
    kk=kk[,id_new]
    K[[i]]=kk
  }
  delta<-datay[,2]
  rho0=rep(-1,N)
  g_rho_l<-rep(0,N)
  g_rho_phi<-rep(0,N)
  g_rho_phim<-list()
  H_rho_phim<-list()
  
  
  for (i in 1:(N-1)) {
    fmu_0=1
    fzi_00=1
    for(j in (i+1):N){
      fmu_1=-sum(rho0[j:N])
      fmu_0=fmu_0*fmu_1
      fzi_1=delta[j]+fmu_1
      fzi_00=fzi_00*fzi_1
    }
    
    g_rho_l[i]=fmu_0/(delta[i]-rho0[i])/fzi_00
    phi_m=t(rho0)%*%K[[i]]%*%rho0
    if(phi_m<=C*(1-lambda)) {g_rho_phim[[i]]=rep(0,N)
    H_rho_phim[[i]]=matrix(0,N,N)
    
    }
    else{
      phi_c=(phi_m-C*(1-lambda))/C/lambda
      g_rho_phim[[i]]=(phi_c/phi_m)*K[[i]]%*%rho0
      u=K[[i]]%*%rho0%*%t(rho0)%*%K[[i]]
      H_rho_phim[[i]]=phi_c*(K[[i]]/phi_m-u/(phi_m^3))+u/(C*lambda)/(phi_m^2)
    }
    
  }
  g_rho_l[N]=-log(delta[N]-rho0[N])
  phi_m=t(rho0)%*%K[[N]]%*%rho0
  if(phi_m<=C*(1-lambda)) {g_rho_phim[[N]]=rep(0,N)
  H_rho_phim[[N]]=matrix(0,N,N)   
  }
  else{
    phi_c=(phi_m-C*(1-lambda))/C/lambda
    g_rho_phim[[N]]=(phi_c/phi_m)*K[[N]]%*%rho0
    u=K[[N]]%*%rho0%*%t(rho0)%*%K[[N]]
    H_rho_phim[[N]]=phi_c*(K[[N]]/phi_m-u/(phi_m^3))+u/(C*lambda)/(phi_m^2)
    
    
  }
  H_rho=matrix(0,N,N)
  for (i in 1:N) {
    g_rho_phi=g_rho_phi+g_rho_phim[[i]]
    H_rho=H_rho+H_rho_phim[[i]]
  }
  g_rho=g_rho_l+g_rho_phi
  hl=1/(delt-rho0)
  H_rho_l=diag(hl)
  H_rho=H_rho+H_rho_l
  
  Am<-matrix(0,n,n)
  for (i in 1:(n-1)) {
    for(j in (i+1):n)
      Am[i,j]=1
    
  }
  Am[n,n]=1
  hg=-solve(H_rho)%*%g_rho
  as1<-Am%*%rho0
  bs1<-Am%*%hg
  cs1<-rep(0,n)
  
  s1<-find_sz(as1,bs1,cs1)
  
  s2<-find_sz(rho0,hg,delta)
  es<-min(s1,s2)
  
  rho1<-rho0+es*hg
}

fit_mk_cox<-function(K,datax,datay,gp,gn,C,lambda,step)
{
  nd<-cox_arr(datax,datay)
  N=length(datax[,1])
  dg<-da_gp(nd[[1]],gn,gp)
  id_new<-nd[[3]]
  #kernel adjust with re-ordered data
  for (m in 1:gn) {
    kk=K[[m]]
    kk=kk[id_new,]
    kk=kk[,id_new]
    K[[m]]=kk
  }
  delta<-datay[,2]
  rho0=rep(-0.1,N)
  g_rho_l<-rep(0,N)
  g_rho_phi<-rep(0,N)
  g_rho_phim<-list()
  H_rho_phim<-list()
  
  rho_up<-function(rho0){
    for (i in 1:(N-1)) {
      fmu_0=1
      fzi_00=1
      for(j in (i+1):N){
        fmu_1=-sum(rho0[j:N])
        fmu_0=fmu_0*fmu_1
        fzi_1=delta[j]+fmu_1
        fzi_00=fzi_00*fzi_1
      }
      g_rho_l[i]=fmu_0/(delta[i]-rho0[i])/fzi_00
    }
    g_rho_l[N]=-log(delta[N]-rho0[N])
    for (m in 1:gn) {
      phi_m=t(rho0)%*%K[[m]]%*%rho0
      phi_m=as.numeric(phi_m)
      print(phi_m)
      if(phi_m<=C*(1-lambda)) {g_rho_phim[[m]]=rep(0,N)
      H_rho_phim[[m]]=matrix(0,N,N)   
      }
      else{
        phi_c=(phi_m-C*(1-lambda))/C/lambda
        g_rho_phim[[m]]=(phi_c/phi_m)*K[[m]]%*%rho0
        u=K[[m]]%*%rho0%*%t(rho0)%*%K[[m]]
        H_rho_phim[[m]]=phi_c*(K[[m]]/phi_m-u/(phi_m^3))+u/(C*lambda)/(phi_m^2)
      }
      
    }
    
    H_rho=matrix(0,N,N)
    for (m in 1:gn) {
      g_rho_phi=g_rho_phi+g_rho_phim[[m]]
      H_rho=H_rho+H_rho_phim[[m]]
    }
    g_rho=g_rho_l+g_rho_phi
    hl=1/(delta-rho0)
    H_rho_l=diag(hl)
    H_rho=H_rho+H_rho_l
    
    Am<-matrix(0,N,N)
    for (i in 1:(N-1)) {
      for(j in (i+1):N)
        Am[i,j]=1
      
    }
    Am[N,N]=1
    hg=-solve(H_rho)%*%g_rho
    as1<-Am%*%rho0
    bs1<-Am%*%hg
    cs1<-rep(0,N)
    
    s1<-find_sz(as1,bs1,cs1)
    print(s1)
    s2<-find_sz(rho0,hg,delta)
    print(s2)
    es<-min(s1,s2)
    
    print(es)
    
    rho1<-rho0+es*hg
    rho1
    print(rho1)
  }
  st=1
  while (st<=step) {
    print(length(rho0))
    rho1<-rho_up(rho0)
    rho0=rho1
    st=st+1
  }
  
  alpha=list()
  w=rep(0,gn)
  for(i in 1:gn){
    phi_m=t(rho0)%*%K[[i]]%*%rho0
    phi_m=as.numeric(phi_m)
    if(phi_m<=C*(1-lambda))
    {
      alpha[[i]]=rep(0,N)
      w[i]=0
    }
    else{
      phi_c=(phi_m-C*(1-lambda))/C/labmdaAllGenes_GenomicRanges.RDS
      WTA_projected.rds
      alpha[[i]]=K[[i]]%*%rho0*phi_c/phi_m
      w[i]=phic/phi_m
    }
    
  }
  
  lis=list(alp=alpha,weight=w,newX=nd[[1]],newY=nd[[2]] )
  
  
  lis
  
  
}