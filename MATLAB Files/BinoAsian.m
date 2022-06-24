%function BinoAsian(S0,E,T,r,sigma,N,F)
%% Test
clear; tic; S0=50;  E=40; T=1.0; r=0.1; sigma=0.3; N=500; F=@(S,A)max(A-E,0);
%% Function to evaluate European Call option by Binomial Method
%   Parameters:
%   S0 = initial share price
%   E = exercise price`
%   T = time to expiry
%   r = riskfree interest rate
%   sigma = volatility
%   N = Number of steps
%   F = Option Payoff (European Call in given)
%% Calculated parameters
dt   = T/N;                  %Timestep
u    = exp(sigma*sqrt(dt));  %Up price movement
d    = 1/u;                  %Down price movement
disf = exp(-r*dt);           %Discount factor over each timestep
p    = (1/disf - d)/(u-d);   %Risk-neutral probability
%% Initalizing Arrays and Functions
S     = zeros(N+1,N+1);   %Underlying Asset Price
S_k   = cell(N+1,N+1);    %S_max %% Cells are used so different sized vectors can be stored at each element in cell array
A     = zeros(N+1,N+1);   %Average of Underlying Asset Price
A_k   = cell(N+1,N+1);    %Representitive averages
C     = zeros(N+1,N+1);   %Price of Option
C_k   = cell(N+1,N+1);    %Option price for given representitive avg
%F     = @(S,A)max(A-E,0); %Option Payoff (European Call)
%% Functions to Calculate Maximum and Minimum Representative Averages
A_min = @(i,j) (1/(i+1)) * (sum(S0*d.^(0:i-j)) + sum(S0*d.^(0:(j-1)+i-(2*j))));
A_max = @(i,j) (1/(i+1)) * (sum(S0*u.^(0:j)) + sum(S0*u.^((0:(i-j-1))+(2*j)-i)));
%% Calculate Underlying asset price
for i = 0:N
    for j = 0:i
        S(i+1,j+1) = S0*u^(j)*d^(i-j);
    end
end
%% Calculating All S_max and Representitive Averages
for i = 0:N
    for j = 0:i %j indexes at j+1 due to matlab not allowing C_k{:,0)
        A_k{i+1,j+1} = zeros(j*(i-j)+1,1); %Create Vector to hold rep avgs
        S_k{i+1,j+1} = NaN(j*(i-j),1);     %Create Vector to hold S_max
        A_k{i+1,j+1}(1) = A_max(i,j);      %Assign A_max to first element in vector
        % Paths with only up (i = j) or down movements (j = 0) or i = 1 will only have one representative average
        if i < 1 || i == j || j == 0 || j*(i-j) == 1
            S_k{i+1,j+1}(1) = S0*(u^j)*d^(i-j);
            continue
        end
        %% Filling S_max Vector
        % Setting Paramaters
        n = j*(i-j);        %Size of vector
        rep = min([j;i-j]); %Maximum element repetition in vector
        % Calculating unique elements in vector
        unique = 2*numel(1:rep-1) + (n - (2*sum(1:rep-1)))/rep;
        % First and last values
        S_k{i+1,j+1}(1)   = S0 * (u ^ j);
        S_k{i+1,j+1}(end) = S0 * (u ^ (j - unique + 1));
        %%The following if statement changes the final output%%
        if n ~= 2
            % Ascending/Descending repeated values
            k     = 2; %iterator
            count = 2; %Number of elements filled
            while k < rep
                S_k{i+1,j+1}(sum(1:k-1)+1:sum(1:k-1)+k) = repmat(S0*(u^(j-k+1)), k, 1);
                S_k{i+1,j+1}(n-sum(1:k-1)+1-k:n-sum(1:k-1)) = repmat(S0*(u^(j-unique+k)), k, 1);
                count = count + 2*k;
                k = k + 1;
            end
            % "Middle" repeated values
            for l = 0 : (n-count)/rep - 1
                S_k{i+1,j+1}(count/2+1+l*rep:end) = [repmat(S0*u^(j-k+1), rep, 1);
                    S_k{i+1,j+1}(count/2+1+l*rep+rep:end)];
                k = k + 1;
            end
        end
        %% Calculating Representitive Averages Using S_max
        A_k{i+1,j+1}(j*(i-j) + 1) = A_min(i,j); %Assign A_min to last element in vector
        for k = 2:j*(i-j)
            A_k{i+1,j+1}(k) = A_k{i+1,j+1}(k-1) - (1/(i+1))* ...
                (S_k{i+1,j+1}(k-1) - S_k{i+1,j+1}(k-1)*(d^2));
        end
    end
end
%% Pricing Option Value at Final Time (N)
for j = 0:N
    C_k{N+1,j+1} = F(S(N+1,j+1),A_k{N+1,j+1});
end
%% Pricing Option
err = 1e-3;
for i = N-1:-1:0
    for j = 0:i
        C_k{i+1,j+1} = zeros(j*(i-j)+1,1);
        for k = 1:j*(i-j)+1
            %% Find K_u
            Ku  = ( (i+1)*A_k{i+1,j+1}(k) + u*S(i+1,j+1) )/(i+2);
            [loc, ubound, lbound] = findInSorted(Ku,A_k{i+2,j+2},err);
            % If found set accordingly
            if loc > 0
                Cu = C_k{i+2,j+2}(loc);
            else
                % If not found, interpolate between closest values
                Cu = C_k{i+2,j+2}(lbound)+(Ku-A_k{i+2,j+2}(lbound))*(...
                    (C_k{i+2,j+2}(ubound)-C_k{i+2,j+2}(lbound))/...
                    (A_k{i+2,j+2}(ubound)-A_k{i+2,j+2}(lbound)));
            end
            %% Find K_d
            Kd = ( (i+1)*A_k{i+1,j+1}(k) + d*S(i+1,j+1) )/(i+2);
            [loc, ubound, lbound] = findInSorted(Kd,A_k{i+2,j+1},err);
            %If found set accordingly
            if loc > 0
                Cd = C_k{i+2,j+1}(loc);
            else
                %If not, interpolate between closest values
                Cd = C_k{i+2,j+1}(lbound)+(Kd-A_k{i+2,j+1}(lbound))*(...
                    (C_k{i+2,j+1}(ubound)-C_k{i+2,j+1}(lbound))/...
                    (A_k{i+2,j+1}(ubound)-A_k{i+2,j+1}(lbound)));
            end
            %% Calculate option value at previous node (i,j decreasing)
            C_k{i+1,j+1}(k) = disf*(p*Cu + (1-p)*Cd);
        end
    end
end
disp(C_k{1,1})
toc;