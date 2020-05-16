// plotCpidResults.sce

M  = fscanfMat("data.txt");
t  = M(:,1);
SP = M(:,2);
PV = M(:,3);
CS = M(:,4);
TS = M(:,5);
TR = M(:,6);
CSo= M(:,7);

scf(0); clf;
subplot(211);plot(t,SP,'k:',t,PV,'r');
subplot(212);plot(t,TR,'k:',t,CS,'b',t,CSo,'b:');

