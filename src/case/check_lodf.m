clc
clear

mpc = ACTIVSg2000;
intmpc = ext2int(mpc);
S = makePTDF(intmpc);
L = abs(makeLODF(intmpc.branch,S));
L(isinf(L)|isnan(L))=0;
L(L<0.01)=0;
L = L-eye(3206,3206);
% imagesc(L)
check1 = max(L,[],1);
rate = intmpc.branch(:,6)/intmpc.baseMVA;

plot(check1,rate,'b*')
hold on

cont = [85; 87; 464; 3007];
crit = check1(cont);
crit_rate = rate(cont);

plot(crit,crit_rate,'r*')
xlabel("Maximum of entries in the LODF corresponding to line")
ylabel("Rating of line")