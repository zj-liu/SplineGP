function [bfm,likv]  = lr(x,data)
%Gaussian spline likelihood ratio matrix
%
% [bfm,likv]  = lr(x,data)
% x: time points column vector
% data: data matrix. each row is the measurement of a gene
%
% bfm: Matrix bfm(i,j) = p(yi,yj)/(p(yi)p(yj))

nRep = 3;
nTime = length(x);
m = length(data);

v = 10.^[-2,-1,0,1,2,3,4,5];
for i = 1:m
    yt = reshape(data(i,:),nRep,nTime);
    yt = yt';
    yt = yt(:);
    for j = 1:length(v)
        likv(i,j) = logml(x,yt,1,nRep,v(j));    
    end
end

v = 100;

for i = 1:m
   yt = reshape(data(i,:),nRep,nTime);
   yt = yt';
   yt = yt(:);
   liki(i) = logml(x,yt,1,nRep,v);
end
likm = zeros(m);
for i = 1:m-1
    parfor j = i+1:m
        yt = reshape(data(i,:),nRep,nTime);
        yt = yt';
        yt = yt(:);
        yt2 = reshape(data(j,:),nRep,nTime);
        yt2 = yt2';
        yt2 = yt2(:);
        yt = [yt;yt2];
        likm(i,j) = logml(x,yt,2,nRep,v);
    end
%     fprintf('gene%d\n',i)
end
likm = likm + likm';
for i = 1:m
    for j = 1:m
        bfm(i,j) = likm(i,j)/liki(i)/liki(j);
    end
end
