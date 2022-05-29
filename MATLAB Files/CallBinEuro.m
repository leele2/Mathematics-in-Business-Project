% function [X, call] = CallBinEuro(S,E,T,N,r,sigma)
clear; tic; S=50;  E=40; T=1; r=0.1; sigma=0.3; N=50;
%CallBinomial.m
%
% Function to evaluate European Call option by Binomial Method
% 
% Parameters are
%    S=initial share price
%    E=exercise price
%    T=time to expiry
%    r=riskfree interest rate
%    sigma=volatility
% Number of steps fixed at 50
% e.g. suitable parameters are  %%%%%
%  S=3;  E=2; T=1; r=0.05; sigma=0.3;

n=N;
% Calculated parameters 
dt=T/n;                 %% Timestep
u=exp(sigma*sqrt(dt));  %% Up price movement
d=1/u;                  %% Down price movement
disf=exp(-r*dt);        %% discount factor over each timestep 
p=(1/disf - d)/(u-d);   %% risk-neutral probability
%%%

X=zeros(n+1,n+1); % Array for option values
                  % Note: MATLAB subscripts start from 1 (not 0)
for k=1:(n+1) 
    X(n+1,k)=max(S*u^(k-1)*d^(n-k+1)-E,0);
end
for m=n:-1:1
    for k=1:m 
        X(m,k)=disf*( p*X(m+1,k+1) + (1-p)*X(m+1,k) );
    end
end

call=X(1,1);