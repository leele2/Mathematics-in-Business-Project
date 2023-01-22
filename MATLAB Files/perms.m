S0 = 10;
u  = 1.1;
d  = 1/u;
N  = 10;
r  = 0.05;
delta_t = 0.1;
p_star = (exp(r*delta_t) - d)/(u-d);
node = zeros(N+1,N+1);
avgs = cell(N+1,N+1);
perm = [];
for i = 0:N
    for j = 0:i
        node(i+1,j+1) = S0*u^j*d^(i-j);
        perm = uniqueperms([ones(1,j), zeros(1,i-j)]);
        if numel(perm) > 0
            for k = 1:numel(perm(:,1))
            avg    = ones(1,numel(perm(k,:))+1);
            avg(1) = S0;
            for l = 1:numel(perm(k,:))
                if perm(k,l) == 1
                    avg(l+1) = avg(l)*u;
                else
                    avg(l+1) = avg(l)*d;
                end
            end
            avgs{i+1,j+1}(k) = mean(avg);
            end
        end
        
    end
end