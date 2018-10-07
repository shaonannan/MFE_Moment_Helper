function [f,g] = optfun(mu,x,lambda0,PEq,fin,gamma)
  % optfun(MU,X,LAMBDA0)
 % This program calculates the Lagrange Multipliers of the ME
 % probability density functions p(x) from the knowledge of the
 % N moment contstraints in the form:
% E{x^n}=mu(n) n=0:N with mu(0)=1.
 %
 % MU is a table containing the constraints MU(n),n=1:N.
 % X is a table defining the range of the variation of x.
% LAMBDA0 is a table containing the first estimate of the LAMBDAs.
 % (This argument is optional.)
 % LAMBDA is a table containing the resulting Lagrange parameters.
 % P is a table containing the resulting pdf p(x).
 % ENTR is a table containing the entropy values at each
 % iteration.
 % f is the value of the function(lambda)
 % g is the gradiant of f
 
 mu=mu(:); mu=[1;mu]; % add mu(0)=1
 x=x(:);  % x axis
 dx=x(2)-x(1);
 lambda=lambda0(:);
 N=length(lambda);
%  fin=zeros(length(x),N); %
%  fin(:,1)=ones(size(x)); % fi0(x)=1
% 
%  for n=2:N
%       fin(:,n)=x.*fin(:,n-1);
%  end
%  p = power(PEq,gamma).*exp((fin(:,1:N)*lambda));  
 p = PEq.*exp((fin(:,1:N)*lambda));  
 
 f = dx*sum(p) - mu'*lambda;
  
 if nargout>1
 G=zeros(N,1); % Calculate Gn
 for n=1:N
  G(n) = trapz(x, fin(:,n).*p);
 end
 g=G(1:N) - mu; 
 end