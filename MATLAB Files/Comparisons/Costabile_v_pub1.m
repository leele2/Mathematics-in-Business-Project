addpath('Functions')  
clear; tic; S0=50;  E=40; T=1; r=0.1; sigma=0.3; N=50; F=@(E,A)max(A-E,0);
F1 = @(E,A) max(mean(A,2)-E,0); M = 1e5;
Ns = [10,15,20:10:90];
BAs = zeros(1,numel(Ns));
BAsP = BAs;
MntC = BAs;
for i = 1:numel(Ns)
    N = Ns(i);
    BAs(i)  = BinoAsian(S0,E,T,r,sigma,N,F);
    MntC(i) = MC(S0,E,T,r,sigma,N,F1,M);
end
pub = [11.5276,11.5348,11.5384,11.5413,11.5302,11.5449,11.5458,11.5463,11.5467,11.547];
vars = {'Alternative Costabile', 'Published Values', 'Monte-Carlo'};
rows = cell(1,numel(Ns));
for i = 1:numel(Ns)
    rows{i} = join(['N = ', int2str(Ns(i))]);
end

table = array2table(round([BAs',pub',MntC'],4),RowNames = rows,VariableNames = vars);
writetable(table,"C:\Users\dj-lu\OneDrive - University of Exeter\University of Exeter\05 - Fifth Year\Mathematics in Business Project\Latex_Files\Main\table2.csv",'WriteRowNames',true)
