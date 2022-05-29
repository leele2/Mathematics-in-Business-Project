tic;
for i = 0:50
    for j = 0:i
        if i < 1 || i == j || j == 0
            out = S0*(u^j)*d^(i-j);
        else
            out = NaN(j*(i-j),1);
            out(1) = S0*(u^j); %S_max for A_max
            %Creates initial path with j up steps followed by i-j downsteps
            %begins at node 1, +x: x steps above s0, -x: x steps below s0
            Tau = [1:j, j-1:-1:j-(i-j)+1];
            %Creates a "minimum" path with i-j downsteps followed by j upsteps
            Tau_min = [-1:-1:-(i-j),-(i-j):-(i-j)+(j-1)];
            %Removes duplicate from min path generation method
            Tau_min(numel(-1:-1:-(i-j))+1) = [];
            for k = 2:j*(i-j)
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
                out(k) = S0 * u^(Tau(Tau_match-1+Index)); %Next S_max is the new Max
            end
        end
    end
end
toc;