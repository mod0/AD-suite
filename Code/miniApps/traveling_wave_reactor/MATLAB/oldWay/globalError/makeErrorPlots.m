function makeErrorPlots(simDir,refDir,docompare,doplots)

load([simDir '/phiStore']);
load([simDir '/timeStore']);
load([simDir '/nineStore']);
load([simDir '/eightStore']);
load([simDir '/zeroStore']);
load([simDir '/SaxStore']);
timeSim=timeStore(1:end-1);
time=timeStore(1:end-1)/(60*60*24);
phiSim=phiStore(1:end-1,:);
N9sim=nineStore(1:end-1,:);
N8sim=eightStore(1:end-1,:);
N0sim=zeroStore(1:end-1,:);
saxSim=SaxStore(1:end-1);
clear phiStore timeStore nineStore zeroStore eightStore SaxStore;

if docompare
    load([refDir '/phiStore']);
    load([refDir '/timeStore']);
    load([refDir '/nineStore']);
    load([refDir '/eightStore']);
    load([refDir '/zeroStore']);
    load([refDir '/SaxStore']);
    timeRef=timeStore(1:end-1);
    phiRef=phiStore(1:end-1,:);
    N9ref=nineStore(1:end-1,:);
    N8ref=eightStore(1:end-1,:);
    N0ref=zeroStore(1:end-1,:);
    saxRef=SaxStore(1:end-1);
    clear phiStore timeStore nineStore zeroStore eightStore SaxStore;
    
    [TF1 LOC1]=ismember(timeRef,timeSim);
    [TF2 LOC2]=ismember(timeSim,timeRef);
    tmptimeRef=timeSim(LOC1(TF1));
    timeSim=timeRef(LOC2(TF2));
    timeRef=tmptimeRef;
    N9ref=N9ref(LOC2(TF2),:);
    phiRef=phiRef(LOC2(TF2),:);
    N8ref=N8ref(LOC2(TF2),:);
    N0ref=N0ref(LOC2(TF2),:);
    saxRef=saxRef(LOC2(TF2));
    
    N9sim=N9sim(LOC1(TF1),:);
    phiSim=phiSim(LOC1(TF1),:);
    N8sim=N8sim(LOC1(TF1),:);
    N0sim=N0sim(LOC1(TF1),:);
    saxSim=saxSim(LOC1(TF1));
    time=timeSim/(60*60*24*365);
end
        

[nTS nCell]=size(phiSim);
h=400/nCell;

getMicros;
sqrtNumPlot=4;
B2g=(pi/150)^2;
nVec=zeros(3,nCell);
cellInds=round(linspace(1,round(2*nCell/3),sqrtNumPlot^2));
alphaNum=zeros(nTS,sqrtNumPlot^2);
alphaInfNum=zeros(nTS,sqrtNumPlot^2);
term1=zeros(nTS,sqrtNumPlot^2);
term2=zeros(nTS,sqrtNumPlot^2);
for ts=1:nTS
    nVec(1,:)=N9sim(ts,:);
    nVec(2,:)=N8sim(ts,:);
    nVec(3,:)=N0sim(ts,:);
    getMacros;
    for cl=1:sqrtNumPlot^2
        thisCell=cellInds(cl);
        if thisCell==1
            alphaNum(ts,cl)=nSpeed*(macro(1,1)-macro(2,1)-...
                saxSim(ts)-macro(3,1)*B2g-...
                bndMeans(1,1)/h^2*(phiSim(ts,1)-phiSim(ts,2))/phiSim(ts,1)-...
                2*macro(3,1)/(h*(h+4*macro(3,1))));
            alphaInfNum(ts,cl)=nSpeed*(macro(1,thisCell)-...
                macro(2,thisCell)-saxSim(ts)-macro(3,thisCell)*B2g);
        elseif thisCell==nCell
            alphaNum(ts,cl)=nSpeed*(macro(1,nCell)-macro(2,nCell)-...
                saxSim(ts)-macro(3,nCell)*B2g-...
                bndMeans(1,nCell-1)/h^2*(phiSim(ts,nCell)-phiSim(ts,nCell-1))/phiSim(ts,nCell)-...
                2*macro(3,nCell)/(h*(h+4*macro(3,nCell))));
            alphaInfNum(ts,cl)=nSpeed*(macro(1,thisCell)-...
                macro(2,thisCell)-saxSim(ts)-macro(3,thisCell)*B2g);
        else
            alphaNum(ts,cl)=nSpeed*(macro(1,thisCell)-...
                macro(2,thisCell)-saxSim(ts)-macro(3,thisCell)*B2g-...
                bndMeans(1,thisCell-1)/h^2*(phiSim(ts,thisCell)-phiSim(ts,thisCell-1))/phiSim(ts,thisCell)-...
                bndMeans(1,thisCell)/h^2*(phiSim(ts,thisCell)-phiSim(ts,thisCell+1))/phiSim(ts,thisCell));
            alphaInfNum(ts,cl)=nSpeed*(macro(1,thisCell)-...
                macro(2,thisCell)-saxSim(ts)-macro(3,thisCell)*B2g);
        end
        term1(ts,cl)=phiSim(ts,thisCell)*...
            (N9sim(ts,thisCell)*(Gamma*msa8+msa9)-N8sim(ts,thisCell)*(msa8+msa8^2/msa9));
        term2(ts,cl)=alphaInfNum(ts,cl)*(msa8/msa9*N8sim(ts,thisCell)-N9sim(ts,thisCell));
    end
end
kill=12;

if doplots
    figure(2)
    clf
    buffer=.025;
    plotDim=(1-buffer*(sqrtNumPlot))/sqrtNumPlot;
    for i=1:sqrtNumPlot
        for j=1:sqrtNumPlot
            thisCell=(i-1)*sqrtNumPlot+j;
            axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlot-i)*(buffer+plotDim) plotDim plotDim])
            plot(time(2:end),(term1(2:end,thisCell)+term2(2:end,thisCell))./(term1(1:end-1,thisCell)+term2(1:end-1,thisCell)))
%             max(abs((term1(1:end-1,thisCell)+term2(1:end-1,thisCell))))
            set(gca,'xtick',[]);set(gca,'ytick',[]);
            set(gca,'xlim',[time(2) time(end)]);
            line([time(2) time(end)],[1 1],'LineStyle','--','Color','k')
            ax=axis;
            tmp=ax(4)-ax(3);
            set(gca,'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
            line([time(2) time(end)],[0 0],'Color','k')
            line([time(2) time(end)],[-1 -1],'LineStyle','--','Color','k')
            set(gca,'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
            text(.1,.9,num2str(cellInds(thisCell)),'units','normalized','FontSize',14)
        end
    end
    
    figure(3)
    clf
    buffer=.025;
    plotDim=(1-buffer*(sqrtNumPlot))/sqrtNumPlot;
    for i=1:sqrtNumPlot
        for j=1:sqrtNumPlot
            thisCell=(i-1)*sqrtNumPlot+j;
            axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlot-i)*(buffer+plotDim) plotDim plotDim])
            plot(time,(N9sim(:,cellInds(thisCell))-N9ref(:,cellInds(thisCell)))./N9ref(:,cellInds(thisCell)))
%             max(abs((N9sim(:,cellInds(thisCell))-N9ref(:,cellInds(thisCell)))./N9ref(:,cellInds(thisCell))))
            set(gca,'xtick',[]);set(gca,'ytick',[]);
            set(gca,'xlim',[time(1) time(end)]);
            ax=axis;
            tmp=ax(4)-ax(3);
            set(gca,'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
            line([time(1) time(end)],[0 0],'LineStyle','--','Color','k')
            text(.1,.9,num2str(cellInds(thisCell)),'units','normalized','FontSize',14)
        end
    end
    
    figure(5)
    clf
    buffer=.025;
    plotDim=(1-buffer*(sqrtNumPlot))/sqrtNumPlot;
    for i=1:sqrtNumPlot
        for j=1:sqrtNumPlot
            thisCell=(i-1)*sqrtNumPlot+j;
            axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlot-i)*(buffer+plotDim) plotDim plotDim])
            [AX H1 H2]=plotyy(time,(N9sim(:,cellInds(thisCell))-N9ref(:,cellInds(thisCell)))./N9ref(:,cellInds(thisCell)),...
                time(2:end),(term1(2:end,thisCell)+term2(2:end,thisCell))./(term1(1:end-1,thisCell)+term2(1:end-1,thisCell)));
            set(AX(1),'YTick',[]); set(AX(1),'XTick',[]);
            set(AX(2),'YTick',[]); set(AX(2),'XTick',[]);
            set(gca,'xlim',[time(1) time(end)]);
            text(.1,.9,num2str(cellInds(thisCell)),'units','normalized','FontSize',14)
            axes(AX(1));
            line([time(1) time(end)],[0 0],'LineStyle','--','Color','b')
            ax=axis;
            tmp=ax(4)-ax(3);
            set(AX(1),'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
            axes(AX(2))
            line([time(2) time(end)],[1 1],'LineStyle','--','Color','g')
            ax=axis;
            tmp=ax(4)-ax(3);
            if (ax(4)<1)
                set(AX(2),'ylim',[ax(3) 1+.1*tmp]);
            elseif ax(3)>1
                set(AX(2),'ylim',[1-.1*tmp ax(4)]);
            end
            set(gca,'xlim',[time(1) time(end)]);
            
%             set(AX(2),'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
%             line([time(2) time(end)],[0 0],'Color','k')
%             line([time(2) time(end)],[-1 -1],'LineStyle','--','Color','k')
%             set(AX(2),'ylim',[ax(3)-.1*tmp ax(4)+.1*tmp]);
        end
    end
    
%     figure(4)
%     clf
%     buffer=.025;
%     plotDim=(1-buffer*(sqrtNumPlot))/sqrtNumPlot;
%     for i=1:sqrtNumPlot
%         for j=1:sqrtNumPlot
%             thisCell=(i-1)*sqrtNumPlot+j;
%             axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlot-i)*(buffer+plotDim) plotDim plotDim])
%             plot(time,term1(:,thisCell)+term2(:,thisCell))
%             set(gca,'xtick',[]);set(gca,'ytick',[]);
%             set(gca,'xlim',[time(1) time(end)]);
%             text(.1,.9,num2str(cellInds(thisCell)),'units','normalized','FontSize',14)
%             line([time(1) time(end)],[term1(1,thisCell)+term2(1,thisCell) term1(1,thisCell)+term2(1,thisCell)],'LineStyle','--')
%             text(time(2),term1(1,thisCell)+term2(1,thisCell),num2str(term1(1,thisCell)+term2(1,thisCell),'%.2e'))
%         end
%     end
end