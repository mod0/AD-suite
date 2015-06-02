close('all')
[nts nCell]=size(phiStore);
xvec=linspace(2,398,nCell);
posInds=find(SaxStore>0);
tyr=timeStore/60/60/24/365;
range=1:10:length(posInds);

figure
axes('FontSize',14)
surf(xvec,tyr(posInds(range)),phiStore(posInds(range),:))
ylabel('Time (yr)')
xlabel('Space (X)')
zlabel('Flux')

figure
axes('FontSize',14)
surf(xvec,tyr(posInds(range)),nineStore(posInds(range),:))
zlabel('Pu-239 Conc')
ylabel('Time (yr)')
xlabel('Space (X)')
set(gca,'Zscale','log')

figure
axes('FontSize',14)
surf(xvec,tyr(posInds(range)),eightStore(posInds(range),:))
ylabel('Time (yr)')
xlabel('Space (X)')
zlabel('U-238 Conc')
set(gca,'Zscale','log')

figure
axes('FontSize',14)
plot(tyr,SaxStore)
xlabel('Time (yr)')
ylabel('External Absorber')
ax=axis;
line([tyr(posInds(end)) tyr(posInds(end))],[ax(3) SaxStore(posInds(end))],...
    'LineStyle','--');
