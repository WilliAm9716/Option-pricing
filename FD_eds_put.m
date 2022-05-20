% Finite Difference - explicit scheme(S) for vanilla put options
%
%   calling syntax:
%   v=FD_eds_put(S0, X, r, T, sig, N, I)
%
function fd_v=FD_eds_put(S0, X, r, T, sig, N, I)
% For Explicit schemes, N  has to be chosen large enough to avoid 
% violating the monotinicity condition
%
Smax=3*X;   % set maximum S to be three times the strike value
dt=T/N;
h=Smax/I;
%
VGrid=zeros(I+1,N+1); % finite difference grid
ishift=1; 

% Boundary conditions
VGrid(1,:)=X*exp(-r*(T-(0:dt:T))); % at S=0;
VGrid(I+1,:)=0; % at Smax, unlikely to have positive payoff for large S

% Terminal condition
VGrid(:,N+1)=max(X-(0:I)*h,0);

i=(1:I-1)';  isq=i.^2;
% Explicit Scheme I
c=(0.5*sig^2*isq+0.5*r*i)*dt;
b=1-sig^2*isq*dt-r*dt;
a=(0.5*sig^2*isq-0.5*r*i)*dt;
% Explicit Scheme II
%  c=(0.5*sig^2*isq+0.5*r*i)*dt/(1+r*dt);
%  b=(1-sig^2*isq*dt)/(1+r*dt);
%  a=(0.5*sig^2*isq-0.5*r*i)*dt/(1+r*dt);
% Explicit Scheme III
% c=(0.5*sig^2*isq+r*i)*dt/(1+r*dt);
% b=(1-sig^2*isq*dt-r*i*dt)/(1+r*dt);
% a=(0.5*sig^2*isq)*dt/(1+r*dt);
% Check on monotonicity
len01=length(find(a<0));
disp(['Coeff a, Of ',num2str(I-1), ' elements, ', num2str(len01),' violated the positivity condition.']);
len02=length(find(b<0));
disp(['Coeff b, Of ',num2str(I-1), ' elements, ', num2str(len02),' violated the positivity condition.']);
len03=length(find(c<0));
disp(['Coeff c, Of ',num2str(I-1), ' elements, ', num2str(len03),' violated the positivity condition.']);

ishift=1;
for n=N:-1:1  % backward time recursive
    VGrid(i+ishift,n)=a.*VGrid(i-1+ishift,n+1)+b.*VGrid(i+ishift,n+1)+c.*VGrid(i+1+ishift,n+1);
end;      
fd_v=VGrid(round(S0/h)+ishift,1);  % nearest point interpolation

ExactValue=Pe(S0,X,r,T,sig,0);
disp(['At S0=',num2str(S0),' exact value=',num2str(ExactValue),' FD value=',num2str(fd_v)]);

function y=Pe(S,X,r,t,sigma,q)
d1=(log(S/X)+(r-q+sigma*sigma/2)*t)/sigma/sqrt(t);
d2=d1-sigma*sqrt(t);
y=X*exp(-r*t)*normcdf(-d2)-S*exp(-q*t)*normcdf(-d1);
   
