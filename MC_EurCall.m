function Euro_Call=MC_EurCall(S0,X,r,T,sigma,q,P)
mu=r-q-sigma^2/2;
epsv=randn(P,1);%vector of P random standard normal
ST=S0*exp(mu*T+epsv*sigma*sqrt(T));
Euro_Call=exp(-r*T)*mean(max(ST-X,0));
