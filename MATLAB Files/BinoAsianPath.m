function value = BinoAsianPath(S0,E,T,r,sigma,N,F)
addpath('Functions')
%% Function to evaluate European Call option by Binomial Method
%   Parameters:
%   S0 = initial share price
%   E = exercise price`
%   T = time to expiry
%   r = riskfree interest rate
%   sigma = volatility
%   N = Number of steps
%   F(E,S) = Option Payoff
%% Calculated parameters
dt   = T/N;                  %% Timestep
u    = exp(sigma*sqrt(dt));  %% Up price movement
d    = 1/u;                  %% Down price movement
disf = exp(-r*dt);           %% Discount factor over each timestep
p    = (1/disf - d)/(u-d);   %% Risk-neutral probability
%% Initalizing Arrays and Functions
S     = zeros(N+1,N+1);   %% Underlying Asset Price
S_k   = cell(N+1,N+1);    %% S_max %% Cells are used so different sized vectors can be stored at each element in cell array
A_k   = cell(N+1,N+1);    %% Representitive averages
C_k   = cell(N+1,N+1);    %% Option price for given representitive avg
%% Functions to Calculate Maximum and Minimum Representative Averages
A_min = @(i,j) (1/(i+1)) * (sum(S0*d.^(0:i-j)) + sum(S0*d.^(0:(j-1)+i-(2*j))));
A_max = @(i,j) (1/(i+1)) * (sum(S0*u.^(0:j)) + sum(S0*u.^((0:(i-j-1))+(2*j)-i)));
%% Calculate Underlying asset price
for i = 0:N
    for j = 0:i
        S(i+1,j+1) = S0*u^(j)*d^(i-j);
    end
end
%% Calculating All Representitive Averages
for i = 0:N
    for j = 0:i %j indexes at j+1 due to matlab not allowing C_k{:,0)
        A_k{i+1,j+1} = zeros(j*(i-j)+1,1); % Create Vector to hold rep avgs
        S_k{i+1,j+1} = zeros(j*(i-j),1); % Create Vector to hold S_max
        A_k{i+1,j+1}(1) = A_max(i,j); %Assign A_max to first element in vector
        %Paths with only up (i = j) or down movements (j = 0) or i = 1 will only have one representative average
        if i == j || j == 0
            S_k{i+1,j+1}(1) = S0*(u^j)*d^(i-j);
            continue
        end
        S_k{i+1,j+1}(1) = S0*(u^j); %S_max for A_max
        %Creates initial path with j up steps followed by i-j downsteps
        %begins at node 1, +x: x steps above s0, -x: x steps below s0
        Tau = [1:j, j-1:-1:j-(i-j)+1];
        %Creates a "minimum" path with i-j downsteps followed by j upsteps
        Tau_min = [-1:-1:-(i-j),-(i-j):-(i-j)+(j-1)];
        %Removes duplicate from min path generation method
        Tau_min(numel(-1:-1:-(i-j))+1) = [];
        A_k{i+1,j+1}(j*(i-j) + 1) = A_min(i,j); %Assign A_min to last element in vector
        for k = 2:j*(i-j)
            A_k{i+1,j+1}(k) = A_k{i+1,j+1}(k-1) - (1/(i+1))* ...
                (S_k{i+1,j+1}(k-1) - S_k{i+1,j+1}(k-1)*(d^2));
            %filter where current path not equal to minimum path
            Tau_filtered = Tau ~= Tau_min;
            %Find first index which is not minimum
            [~, Tau_match] = find(Tau_filtered,1);
            %Index highest point on path that isn't on minimum path
            [~, Index] = max(Tau(Tau_filtered));
            %minus two from highest point not on minimum path
            Tau(Tau_match-1+Index) = Tau(Tau_match-1+Index) - 2;
            %Find new highest point not on min path by repeating above steps
            %1: Filter from not being on minimum path
            Tau_filtered = Tau ~= Tau_min;
            %2: Find first index not on minimum path
            [~, Tau_match] = find(Tau_filtered,1);
            %3: Find maximum not on minimum path
            [~, Index] = max(Tau(Tau_filtered));
            %Use new highest point to calculate the next S_max value
            S_k{i+1,j+1}(k) = S0 * u^(Tau(Tau_match-1+Index)); %Next S_max is the new Max
        end
    end
end
%% Pricing Option Value at Final Time (N)
for j = 0:N
    C_k{N+1,j+1} = F(E,A_k{N+1,j+1});
end
%% Pricing Option
err = 1e-6;
for i = N-1:-1:0
    for j = 0:i
        C_k{i+1,j+1} = zeros(j*(i-j)+1,1);
        for k = 1:j*(i-j)+1
            %Find K_u
            Ku  = ( (i+1)*A_k{i+1,j+1}(k) + u*S(i+1,j+1) )/(i+2);
            [loc, ubound, lbound] = findInSorted(Ku,A_k{i+2,j+2},err);
            if loc > 0
                % Set value if found
                Cu = C_k{i+2,j+2}(loc);
            else
                % Otherwise, interpolate
                Cu = C_k{i+2,j+2}(lbound)+(Ku-A_k{i+2,j+2}(lbound))*(...
                    (C_k{i+2,j+2}(ubound)-C_k{i+2,j+2}(lbound))/...
                    (A_k{i+2,j+2}(ubound)-A_k{i+2,j+2}(lbound)));
            end
            %Find K_d
            Kd = ( (i+1)*A_k{i+1,j+1}(k) + d*S(i+1,j+1) )/(i+2);
            [loc, ubound, lbound] = findInSorted(Kd,A_k{i+2,j+1},err);
            if loc > 0
                % Set value if found
                Cd = C_k{i+2,j+1}(loc);
            else
                % Otherwise, interpolate
                Cd = C_k{i+2,j+1}(lbound)+(Kd-A_k{i+2,j+1}(lbound))*(...
                    (C_k{i+2,j+1}(ubound)-C_k{i+2,j+1}(lbound))/...
                    (A_k{i+2,j+1}(ubound)-A_k{i+2,j+1}(lbound)));
            end
            % Price option using discounted expected value
            C_k{i+1,j+1}(k) = disf*(p*Cu + (1-p)*Cd);
        end
    end
end
value = C_k{1,1};