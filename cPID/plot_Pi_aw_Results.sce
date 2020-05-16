// plotCpidResults.sce

M  = fscanfMat("datair.txt");
t  = M(:,1);
SP = M(:,2);
PV = M(:,3);
CS = M(:,4);
TS = M(:,5);
TR = M(:,6);

scf(0); clf;
subplot(211);plot(t,SP,'k:',t,PV,'r');
subplot(212);plot(t,TR,'k:',t,CS,'b');

