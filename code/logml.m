function lik = logml(x,y,nGene,nRep,v)
%Likelihood of given genes using the Gaussian process model and 1-order
%truncated power spline basis.
%
% lik = logml(x,y,nGene,nRep,v)
% x: time points column vector
% y: gene data, matrix of size (nRep*Time)*nGene 
%       each column is of length (nRep*Time) and in the form of
%       [y(1,t1),...y(1,tT),y(2,t1),...,y(2,tT),...y(nRep,t1),...y(nRep,tT)]'
% nRep:  number of replicates
% nGene:  number of genes
% v: prior covariance

%%TODO input parameter check
y = y(:);
N = nGene*nRep;

T = length(x);
p = length(x);
q = 1; %spline order

%%%parameter
Alpha = 0.01; %
Gamma = 0.01; %
% v = 1; %

X = zeros(T,T);
V = v*eye(p,p);

g = sqrt(Gamma^Alpha/pi^(N*T))*gamma((N*T+Alpha)/2)/gamma(Alpha/2);

X = zeros(T,T);
for i = 1:T
    X(i,1) = 1;
    X(i,2:end) = max((x(i) - x(1:end-1)'),0).^q;
end
X1 = repmat(X,N,1);

invVp = X1'*X1+inv(V);
m = invVp\X1'*y;
d = y'*y - y'*X1*m;

lik = g/sqrt(det(invVp)) / sqrt(det(V)) / (d+Gamma)^((N*T+Alpha)/2);