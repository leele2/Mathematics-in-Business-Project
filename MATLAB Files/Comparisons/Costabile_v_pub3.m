addpath('Functions')  
clear; tic; S0=100;  E=100; T=0.25; r=0.1; sigma=0.1; N=50; F=@(E,A)max(A-E,0);
F1 = @(E,A) max(mean(A,2)-E,0); M = 1e6;
Ns = [10,15,20:10:90];
BAs = zeros(1,numel(Ns));
BAsP = BAs;
MntC = BAs;
for i = 1:numel(Ns)
    N = Ns(i);
    BAs(i)  = BinoAsian(S0,E,T,r,sigma,N,F);
    MntC(i) = MC(S0,E,T,r,sigma,N,F1,M);
end
pub = [1.8381,1.8418,1.8442,1.8466,1.8475,1.8485,1.849,1.8492,1.8497,1.8499];
vars = {'Alternative Costabile', 'Published Values', 'Monte-Carlo'};
rows = cell(1,numel(Ns));
for i = 1:numel(Ns)
    rows{i} = join(['N = ', int2str(Ns(i))]);
end

table = array2table(round([BAs',pub',MntC'],4),RowNames = rows,VariableNames = vars);
writetable(table,"C:\Users\dj-lu\OneDrive - University of Exeter\University of Exeter\05 - Fifth Year\Mathematics in Business Project\Latex_Files\Main\table4.csv",'WriteRowNames',true)