% Sample BTM program for European vanilla call options
% call syntax: OptVal=btm_EurCall(S0,X,r,T,sigma,q,N)
function OptVal=btm_EurCall(S0,X,r,T,sigma,q,N)
% step 1: set up tree parameters
dt=T/N;                       % delta t, length of time period
u=exp(sigma*sqrt(dt)); d=1/u; % up and down factors
p=(exp((r-q)*dt)-d)/(u-d);    % risk-neutral probability
% step 2: terminal condition 
j = 0:1:N;   % range of index for price states at expiry 
Vn=max(S0*u.^(2*j-N)-X,0);
% step 3: backward recursive through time
jshift = 1;
for n=N-1:-1:0   
   Vnp1=Vn;
   j = 0:1:n; 
   Vn=exp(-r*dt)*(p*Vnp1(j+1+jshift)+(1-p)*Vnp1(j+jshift));
end
% step 4: 
OptVal=Vn(0+jshift);


