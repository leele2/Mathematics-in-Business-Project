N = 50;
test = cell(N+1);
tic;
for i = 0:50
    for j = 0:i
        if i < 1 || i == j || j == 0 || j*(i-j) == 1
            test{i+1,j+1} = S0*(u^j)*d^(i-j);
            if sum(test{i+1,j+1} - S_k{i+1,j+1}) ~= 0
                disp("broke")
            end
        else
            %% Setting Paramaters
            n = j * (i - j);         %Size of vector
            r = min([j; i - j]);     %Maximum element repetition in vector
            out = NaN(n, 1);         %Creating output vector
            count = 2;
            %% Filling S_max vector
            %Calculates unique elements in vector (helps in minimum value calc)
            unique = 2*numel(1:r-1) + (n - (2*sum(1:r-1)))/r;
            %First and last values
            out(1)   = S0 * (u ^ j);
            out(end) = S0 * (u ^ (j - unique + 1));
%             if n == 2
%                 continue
%             end
            %Ascending/Descending repeated values
            k = 2;
            while k < r
                out(sum(1:k-1)+1:sum(1:k-1)+k) = repmat(S0*(u^(j-k+1)), k, 1);
                out(n-sum(1:k-1)+1-k:n-sum(1:k-1)) = repmat(S0*(u^(j-unique+k)), k, 1);
                count = count + 2*k;
                k = k + 1;
            end
            %"Middle" repeated values
            for l = 0 : (n-count)/r - 1
                out(count/2+1+l*r:end) = [repmat(S0*u^(j-k+1), r, 1);
                                          out(count/2+1+l*r+r:end)];
                k = k + 1;
            end
            test{i+1,j+1} = out;
            if sum(test{i+1,j+1} - S_k{i+1,j+1}) ~= 0
                disp("broke")
            end
        end
    end
end
toc;