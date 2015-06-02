function determineOrderInTime(doplots)

rng=1:400;

load('12hr/phiStore');
load('12hr/timeStore');
load('12hr/nineStore');
load('12hr/eightStore');
load('12hr/zeroStore');
load('12hr/SaxStore');
time12=timeStore(1:end-1)/(60*60*24);
phi12=phiStore(1:end-1,:);
N912=nineStore(1:end-1,:);
N812=eightStore(1:end-1,:);
N012=zeroStore(1:end-1,:);
sax12=SaxStore(1:end-1);
clear phiStore timeStore nineStore zeroStore eightStore SaxStore;


load('6hr/phiStore');
load('6hr/timeStore');
load('6hr/nineStore');
time06=timeStore(1:2:end);
phi06=phiStore(1:2:end,:);
N906=nineStore(1:2:end,:);
clear phiStore timeStore nineStore zeroStore eightStore SaxStore;

load('3hr/phiStore');
load('3hr/timeStore');
load('3hr/nineStore');
load('3hr/eightStore');
load('3hr/zeroStore');
load('3hr/SaxStore');
time03=timeStore(1:4:end);
phi03=phiStore(1:4:end,:);
N903=nineStore(1:4:end,:);
N803=eightStore(1:4:end,:);
N003=zeroStore(1:4:end,:);
sax03=SaxStore(1:4:end);
clear phiStore timeStore nineStore zeroStore eightStore SaxStore;

nTS=length(time12);
nCell=size(phi03,2);
h=400/nCell;
phiOrder=zeros(nTS);
nineOrder=phiOrder;
vPhi=zeros(nTS,length(rng));
masterTS=12*60*60;
for i=1:nTS
    tmp=norm(phi12(i,rng)-phi06(i,rng))/norm(phi06(i,rng)-phi03(i,rng));
    phiOrder(i)=log(tmp)/log(2);
    vPhi(i,:)=(phi12(i,rng)-phi03(i,rng))/(masterTS^phiOrder(i)-(masterTS/4)^phiOrder(i));
    tmp=norm(N912(i,rng)-N906(i,rng))/norm(N906(i,rng)-N903(i,rng));
    nineOrder(i)=log(tmp)/log(2);
end

kill=12;

if doplots
    figure(1)
    clf
    axes('FontSize',14)
    plot(time12,phiOrder,'r--',time12,nineOrder,'k^');
    xlabel('Time (d)');
    ylabel('Estimated Order in Time')
    legend('\phi Metric','N_9 Metric')
end




% getMicros;
% sqrtNumPlot=sqrt(10);
% B2g=(pi/150)^2;
% nVec=zeros(3,nCell);
% cellInds=round(linspace(1,nCell,sqrtNumPlot^2));
% cellInds=196:205;
% alphaNum=zeros(nTS,sqrtNumPlot^2);
% term1=zeros(nTS,sqrtNumPlot^2);
% term2=zeros(nTS,sqrtNumPlot^2);
% for ts=1:nTS
%     nVec(1,:)=N912(ts,:);
%     nVec(2,:)=N812(ts,:);
%     nVec(3,:)=N012(ts,:);
%     getMacros;
%     for cl=1:sqrtNumPlot^2
%         thisCell=cellInds(cl);
%         if thisCell==1
%             alphaNum(ts,cl)=nSpeed*(macro(1,1)-macro(2,1)-...
%                 sax12(ts)-macro(3,1)*B2g-...
%                 bndMeans(1,1)/h^2*(phi12(ts,1)-phi12(ts,2))/phi12(ts,1)-...
%                 2*macro(3,1)/(h*(h+4*macro(3,1))));
%         elseif thisCell==nCell
%             alphaNum(ts,cl)=nSpeed*(macro(1,nCell)-macro(2,nCell)-...
%                 sax12(ts)-macro(3,nCell)*B2g-...
%                 bndMeans(1,nCell-1)/h^2*(phi12(ts,nCell)-phi12(ts,nCell-1))/phi12(ts,nCell)-...
%                 2*macro(3,nCell)/(h*(h+4*macro(3,nCell))));
%         else
%             alphaNum(ts,cl)=nSpeed*(macro(1,thisCell)-...
%                 macro(2,thisCell)-sax12(ts)-macro(3,thisCell)*B2g-...
%                 bndMeans(1,thisCell-1)/h^2*(phi12(ts,thisCell)-phi12(ts,thisCell-1))/phi12(ts,thisCell)-...
%                 bndMeans(1,thisCell)/h^2*(phi12(ts,thisCell)-phi12(ts,thisCell+1))/phi12(ts,thisCell));
%         end
%         term1(ts,cl)=phi12(ts,thisCell)*...
%             (N912(ts,thisCell)*(Gamma*msa8+msa9)-N812(ts,thisCell)*(msa8+msa8^2/msa9));
%         term2(ts,cl)=alphaNum(ts,cl)*(msa8/msa9*N812(ts,thisCell)-N912(ts,thisCell));
%     end
% end
% kill=12;
% 
% if doplots
%     figure(2)
%     clf
%     buffer=.025;
%     plotDim=(1-buffer*(sqrtNumPlot))/sqrtNumPlot;
%     for i=1:sqrtNumPlot
%         for j=1:sqrtNumPlot
%             thisCell=(i-1)*sqrtNumPlot+j;
%             axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlot-i)*(buffer+plotDim) plotDim plotDim])
%             plot(time12(2:end),(term1(2:end,thisCell)+term2(2:end,thisCell))./(term1(1:end-1,thisCell)+term2(1:end-1,thisCell)))
%             set(gca,'xtick',[]);set(gca,'ytick',[]);
%             set(gca,'xlim',[time12(2) time12(end)]);
%             line([time12(2) time12(end)],[1 1],'LineStyle','--','Color','k')
% %             line([time12(2) time12(end)],[0 0],'Color','k')
% %             line([time12(2) time12(end)],[-1 -1],'LineStyle','--','Color','k')
% %             [AX H1 H2]=plotyy(time(timeInds),[errOne(:,thisCell) errTwo(:,thisCell)],time(timeInds),alphaEigInf(:,thisCell));
% %             set(AX(1),'YTick',[]); set(AX(1),'XTick',[]);
% %             set(AX(2),'YTick',[]); set(AX(2),'XTick',[]);
% %             set(AX(1),'XLim',[0 max(time)]); set(AX(2),'XLim',[0 max(time)]);
% 
% 
%             text(.1,.9,num2str(cellInds(thisCell)),'units','normalized','FontSize',14)
%         end
%     end
% end




