function [sigmas,Wm,Wc]=SigmaPoints(n,alpha,beta,kappa,x,P)
% Generates sigma points and weights with the given mean and covariance
% x : Given Mean 
% P : Given Covariance 
% n : dimension of the state (2n+1 weights will be generated)
% alpha : the spread of the sigma points around the mean. (usually a small positive value)
% beta : incorporates prior knowlegdge of the distribution of the mean
% kappa : Secondary scaling parameter uduslly set to 0 or to 3-n
% sigmas n*(2n+1) :array of sigma points
% Wm (2n+1)*1: weights of mean
% Wc (2n+1)*1：weights of covariance

% Generate sigma points
lambda= alpha^2*(n+kappa)-n;
U=chol(P*(lambda+n));
sigmas=zeros(2*n+1,n);
sigmas(1,:)=x;
for i=1:n%2:2*n+1 
    sigmas(2*i,:)=x-U(i);
    sigmas(2*i+1,:)=x+U(i);
end
% Compute weights
Wm=zeros(2*n+1,1)+0.5/(n+lambda);
Wc=zeros(2*n+1,1)+0.5/(n+lambda);
Wm(1)=lambda/(n+lambda);
Wc(1)=lambda/(n+lambda)+(1-alpha^2+beta);
end

