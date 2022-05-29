tic;
test = cell(51,51);
for i = 0:50
    for j = 0:i
        if i < 1 || i == j || j == 0 || j*(i-j) == 1
            out = S0*(u^j)*d^(i-j);
        else
            %% Setting Paramaters
            n = j * (i - j);         %Size of vector
            r = min([j; i - j]);     %Maximum element repetition in vector
            out = NaN(n, 1);         %Creating output vector
            %% Filling S_max vector
            %Calculates unique elements in vector (helps in minimum value calc)
            unique = 2*numel(1:r-1) + (n - (2*sum(1:r-1)))/r;
            %First and last values
            out(1)   = S0 * (u ^ j);
            out(end) = S0 * (u ^ (j - unique + 1));
            if n == 2
                continue
            end
            %Ascending/Descending repeated values
            k = 2;
            while k < r
                out(sum(1:k-1)+1:sum(1:k-1)+k) = repmat(S0*(u^(j-k+1)), k, 1);
                out(n-sum(1:k-1)+1-k:n-sum(1:k-1)) = repmat(S0*(u^(j-unique+k)), k, 1);
                k = k + 1;
            end
            %"Middle" repeated values
            c = sum(2:k-1);
            for l = 0 : (n - (2*sum(1:r-1)))/r - 1
                out(2 + c + l*r:end) = [repmat(S0*u^(j-k+1), r, 1); out(2+r+c + l*r:end)];
                k = k + 1;
            end
            if sum(out - S_k{i+1,j+1}) ~= 0
                disp("broke")
            end
        end
        test{i+1,j+1} = out;
    end
end
toc;