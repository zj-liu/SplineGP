function yy = gpsplineregrssion(x,y,nGene,nRep,v)
%Bayesian regression using the Gaussian process model and 1-order
%truncated power spline basis.
%
% yy = gpsplineregrssion(x,y,nGene,nRep,v)
% x: time points column vector
% y: gene data, matrix of size (nRep*Time)*nGene 
%       each column is of length (nRep*Time) and in the form of
%       [y(1,t1),...y(1,tT),y(2,t1),...,y(2,tT),...y(nRep,t1),...y(nRep,tT)]'
% nRep:  number of replicates
% nGene:  number of genes
% v: prior covariance

y = y(:);
%%%parameter
q = 1; %spline order
Alpha = 0.01; %
Gamma = 0.01; %
% v = 100000;

N = nGene*nRep;

T = length(x);
p = length(x);
X = zeros(T,T);
for i = 1:T
    X(i,1) = 1;
    X(i,2:end) = max((x(i) - x(1:end-1)'),0).^q;
end

V = v*eye(p,p);
X1 = repmat(X,N,1);
Vp = inv(X1'*X1+inv(V));
m = Vp*X1'*y;
xx = linspace(min(x),max(x),100)';
XX = zeros(length(xx),T);
for i = 1:length(xx)
    XX(i,1) = 1;
    XX(i,2:end) = max((xx(i) - x(1:end-1)'),0).^q;
end
yy = XX*m; %regression
xp = repmat(x,N,1);

% the variance interval
invVp = X1'*X1+inv(V);
m = invVp\X1'*y;
d = y'*y - y'*X1*m;
sigmastar = (d+Gamma)/(N*T+Alpha-2);
Vy = sigmastar*(XX*Vp*XX'+1);
Vy = diag(Vy);

% plot the interval of +/- 2*std
f = [yy+2*sqrt(Vy(:)); flip(yy-2*sqrt(Vy(:)),1)];
fill([xx; flip(xx,1)], f, [7 7 7]/8);hold on,
for i = 1:N
    plot(x,y(1+(i-1)*length(x):length(x)+(i-1)*length(x)),'*');
end
plot(xx,yy,'r'),ylim([-1,1])
ylim([-max(abs(y)),max(abs(y))])

%%
    function reshapedata(in)
        ids = [1,2];
        ytmp = [];
        for i = 1:length(ids)      
            yt = reshape(data(ids(i),:),3,5);
            yt = yt';
            yt = yt(:);     
            ytmp = [ytmp;yt];
        end   
        figure,gpsplineregrssion(x,ytmp,length(ids),3,100)
    end

end