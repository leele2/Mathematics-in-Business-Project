%function BinoAsian(S0,E,T,r,sigma,N)
%% Test
clear; tic; S0=50;  E=40; T=1; r=0.1; sigma=0.3; N=50;
%% Function to evaluate European Call option by Binomial Method
% Parameters are
%   S0=initial share price
%   E=exercise price
%   T=time to expiry
%   r=riskfree interest rate
%   sigma=volatility
%   N=Number of steps
%% Calculated parameters
dt   = T/N;                  %% Timestep
u    = exp(sigma*sqrt(dt));  %% Up price movement
d    = 1/u;                  %% Down price movement
disf = exp(-r*dt);           %% Discount factor over each timestep
p    = (1/disf - d)/(u-d);   %% Risk-neutral probability
%alpha = 1;                  %% Determines size of representative averages
%h0   = alpha*sigma*sqrt(T); %% Used in caculating h
%h    = h0/(1+(N/100));      %% Used to calculate representitive averages
%% Initalizing Arrays and Functions
S     = zeros(N+1,N+1);  %% Underlying Asset Price
S_k   = cell(N+1,N+1);   %% S_max %% Cells are used so different sized vectors can be stored at each element in cell array
A     = zeros(N+1,N+1);  %% Average of Underlying Asset Price
A_k   = cell(N+1,N+1);   %% Representitive averages
C     = zeros(N+1,N+1);  %% Price of Option
C_k   = cell(N+1,N+1);   %% Option price for given representitive avg
claim = @(x) max(x-E,0); %% Claim function (Call Option)
%% Functions to Calculate Maximum and Minimum Representative Averages
A_min = @(i,j) (1/(i+1)) * (sum(S0.*d.^[0:i-j]) + sum(S0.*d.^([0:(j-1)]+i-(2*j))));
A_max = @(i,j) (1/(i+1)) * (sum(S0.*u.^[0:j]) + sum(S0.*u.^([0:(i-j-1)]+(2*j)-i)));
%% Calculate Underlying asset price
for i = 0:N
    for j = 0:i
        S(i+1,j+1) = S0*u^(j)*d^(i-j);
    end
end
%% Calculating All Representitive Averages
for i = 1:N+1
    for j = 1:i+1 %j indexes at j+1 due to matlab not allowing C_k{:,0)
        %%Debuging
%         i = 10; j = 5;
        A_k{i,j} = zeros((j-1)*(i-(j-1))+1,1); % Create Vector to hold rep avgs
        S_k{i,j} = zeros(((j-1)*(i-(j-1))),1); % Create Vector to hold S_max
        A_k{i,j}(1) = A_max(i,j-1); %Assign A_max to first element in vector
        %Paths with only up (i = j) or down movements (j = 0) or i = 1 will only have one representative average
        if i ~= 1 && i ~= j-1 && j-1 ~= 0
            S_k{i,j}(1) = S0*(u^(j-1)); %S_max for A_max
            Tau = [1:j-1, j-2:-1:(j-1)-(i-(j-1))+1];
            Tau_min = [-1:-1:-(i-(j-1)),-(i-(j-1)):-(i-(j-1))+(j-1)-1];
            Tau_min(numel(-1:-1:-(i-(j-1)))+1) = [];
            %             S_k{i,j}((j-1)*(i-(j-1)) + 1) = S0*u^(2*(j-1)-i); %S_min for A_min %%This is wrong
            A_k{i,j}((j-1)*(i-(j-1)) + 1) = A_min(i,j-1); %Assign A_min to last element in vector
            for k = 2:(j-1)*(i-(j-1))
                A_k{i,j}(k) = A_k{i,j}(k-1) - (1/(i+1))* ...
                    (S_k{i,j}(k-1) - S_k{i,j}(k-1)*(d^2));
                %                 %%Test new method of finding S_max
                %                 [~, Index] = sort(Tau,'descend');
                %                 %Ensuring path doesn't double jump
                %                 for p = 1:length(Tau)
                %                     Tau_test = Tau;
                %                     Tau_test(Index(p)) = Tau(Index(p)) - 2;
                % %                     if abs(max(Tau_test(1:length(Tau)-1) - Tau_test(2:length(Tau)))) > 1
                %                     if abs(max(Tau_test(1:end-1) - Tau_test(2:end))) > 1
                %                         continue
                %                     else
                %                         Tau = Tau_test
                %                         break
                %                     end
                %                     end
                %%New method of finding S_max using Tau_min
%                 if k == 19
%                     disp('here')
%                 end
                Tau_filtered = Tau ~= Tau_min;
                [~, Tau_match] = find(Tau_filtered,1);
%                 Tau_match = numel(Tau) - sum(Tau_filtered);
                [~, Index] = max(Tau(Tau_filtered));
                Tau(Tau_match-1+Index) = Tau(Tau_match-1+Index) - 2;
%                 Tau
                Tau_filtered = Tau ~= Tau_min;
                [~, Tau_match] = find(Tau_filtered,1);
%                 Tau_match = numel(Tau) - sum(Tau_filtered);
                [~, Index] = max(Tau(Tau_filtered));
%                 Tau(Tau_match-1+Index);
                S_k{i,j}(k) = S0 * u^(Tau(Tau_match-1+Index)); %Next S_max is the new Max
            end
        end
    end
end
% %% Pricing Options
% for i = N:-1:0
%     for j = 1:0
%         i = i-1; j = j-1;
%         if j <= i
%             C_k{i,j+1} = zeros(j*(i-j),1);
%             for k = 1:j*(i-j)+1
%                 %Find K_u and K_d
%                 Ku = ( (i+1)*A_k{i,j+1}(k) + u*S(i,j+1) )/(i+2);
%                 if find(A_k{i+1,j+1} == Ku,1) > 0
%                     disp('It was found')
%                 else
%                     Ku = interp1(sort(A_k{i+1,j+1}),Ku,'linear')
%                 end
%             Kd = ( (i+1)*A_k{i,j+1}(k) + d*S(i,j+1) )/(i+2);
%             C_k{i,j+1}(k) = disf*(p*Ku + (1-p)*Kd);
%             end
%         end
%     end
% end
toc;