% Numerical Inversion of Laplace Transforms from INVLAP
% by Juraj 
% https://www.mathworks.com/matlabcentral/fileexchange/32824-numerical-inversion-of-laplace-transforms-in-matlab
function [ft]=inverseLaplace(Fs,radt,a,ns,nd)
% Fs is formula for F(s) as a string
% tini, tend are limits of the solution interval
% nnt is total number of time instants
% a, ns, nd are parameters of the method 
% if not given, the method uses implicit values a=6, ns=20, nd=19
% it is recommended to preserve a=6
% increasing ns and nd leads to lower error
% an example of function calling  
% [t,ft]=INVLAP('s/(s^2+4*pi^2)',0,10,1001);
% to plot the graph of results write plot(t,ft), grid on, zoom on
%FF=strrep(strrep(strrep(Fs,'*','.*'),'/','./'),'^','.^');
if nargin==2
  a=6; ns=20; nd=19;
end    % implicit parameters
%radt=linspace(tini,tend,nnt); % time vector
if radt(1)==0
    radt=radt(2:end);
end  % t=0 is not allowed
alfa = zeros(1,ns+1+nd);
beta = zeros(1,ns+1+nd);
for n=1:ns+1+nd               % prepare necessary coefficients
   alfa(n)=a+(n-1)*pi*1j;
   beta(n)=-exp(a)*(-1)^n;
end
n=1:nd;
bdif=fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
beta(ns+2:ns+1+nd)=beta(ns+2:ns+1+nd).*bdif;
beta(1)=beta(1)/2;
ft = zeros(1,length(radt));
for kt=1:length(radt)                  % cycle for time t
   tt=radt(kt);
   s=alfa/tt;                 % complex frequency s
   bt=beta/tt;
   btF=bt.*feval(Fs,s);          % functional value F(s)
   ft(kt)=sum(real(btF));     % original f(tt)
end
