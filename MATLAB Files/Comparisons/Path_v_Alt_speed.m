addpath('Functions')  
clear; tic; S0=50;  E=40; T=1; r=0.1; sigma=0.3; N=50; F=@(E,A)max(A-E,0);
F1 = @(E,A) max(mean(A,2)-E,0); M = 1e6;
Ns = [10,15,20:10:100];
BAs = zeros(1,numel(Ns));
BAsP = BAs;
for i = 1:numel(Ns)
    N = Ns(i);
    f = @() BinoAsian(S0,E,T,r,sigma,N,F);
    BAs(i)  = timeit(f);
    f = @() BinoAsianPath(S0,E,T,r,sigma,N,F);
    BAsP(i) = timeit(f);
end
vars = {'Alternative Costabile', 'Path Costabile'};
rows = cell(1,numel(Ns));
for i = 1:numel(Ns)
    rows{i} = join(['N = ', int2str(Ns(i))]);
end
table = array2table(round([BAs',BAsP'],4),RowNames = rows,VariableNames = vars);
writetable(table,"C:\Users\dj-lu\OneDrive - University of Exeter\University of Exeter\05 - Fifth Year\Mathematics in Business Project\Latex_Files\Main\table5.csv",'WriteRowNames',true)