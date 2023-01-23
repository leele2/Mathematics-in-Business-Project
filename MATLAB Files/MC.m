function value = MC(S0,E,T,r,sigma,N,F,M)
%% Evaluate Asian type option using Monte-Carlo
%   Parameters:
%   S0 = initial share price
%   E = exercise price`
%   T = time to expiry
%   r = riskfree interest rate
%   sigma = volatility
%   N = number of steps
%   F(E,S) = payoff function
%   M = number of realizations;
%% Calculated parameters
dt   = T/N;       %% Timestep
disf = exp(-r*T); %% Discount factor over each timestep
%% Initalizing Arrays and Functions
S = zeros(M, N+1);
Z = randn(M, N+1);
S(:,1) = S0;
%% Functions to Simulate Realisation
for i = 1:N
    S(:,i+1) = S(:,i).*exp((r-(.5*sigma^2))*dt + (sigma*sqrt(dt).*Z(:,i)));
end
%% Pricing Option
value = disf*mean(F(E, S));