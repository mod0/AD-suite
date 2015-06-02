function compareSpatialDisc


L=400;
h=L/40;
R=150;
B2g=pi^2/R^2;

cd 400cell_3yr
load eightStore
load nineStore
load zeroStore
load SaxStore
load phiStore
load timeStore;
cd ..

ref8=eightStore; ref9=nineStore; ref0=zeroStore; refPhi=phiStore; 
refSax=SaxStore; time=timeStore/(60*60*24*365);
clear eighStore nineStore zeroStore SaxStore phiStore;

cd 40cell_3yr
load eightStore
load nineStore
load zeroStore
load SaxStore
load phiStore
cd ..

one8=eightStore; one9=nineStore; one0=zeroStore; onePhi=phiStore; 
oneSax=SaxStore;
clear eighStore nineStore zeroStore SaxStore phiStore;

cd 40cell_3yr_srtWithRef
load eightStore
load nineStore
load zeroStore
load SaxStore
load phiStore
cd ..

two8=eightStore; two9=nineStore; two0=zeroStore; twoPhi=phiStore; 
twoSax=SaxStore;
clear eighStore nineStore zeroStore SaxStore phiStore;

getMicros;
nVec=zeros(3,40);
nTS=size(ref8,1);
nPlotTimes=100;
timeInds=round(linspace(1,nTS,nPlotTimes));
alphaEigInf=zeros(nPlotTimes,40);
alphaEigNum=zeros(nPlotTimes,40);
ref409=zeros(nPlotTimes,40);
refPhi40=zeros(nPlotTimes,40);
for kk=1:nPlotTimes
    ii=timeInds(kk);
    for j=1:40
        nVec(1,j)=mean(ref9(ii,(j-1)*10+1:j*10));
        nVec(2,j)=mean(ref8(ii,(j-1)*10+1:j*10));
        nVec(3,j)=mean(ref0(ii,(j-1)*10+1:j*10));
        refPhi40(kk,j)=mean(refPhi(ii,(j-1)*10+1:j*10));
    end
    ref409(kk,:)=nVec(1,:);
    getMacros;
    alphaEigInf(kk,:)=macro(1,:)-macro(2,:)-B2g*macro(3,:)-refSax(ii);
    alphaEigNum(kk,2:39)=alphaEigInf(kk,2:39)-...
        bndMeans(1,1:38)/h^2*(1-refPhi40(kk,1:38)/refPhi40(kk,2:39))-...
        bndMeans(1,2:39)/h^2*(1-refPhi40(kk,3:40)/refPhi40(kk,2:39));
    alphaEigNum(kk,1)=alphaEigInf(kk,1)-...
        bndMeans(1,1)/h^2*(1-refPhi40(kk,2)/refPhi40(kk,1))-...
        2*macro(3,1)/(h*(h+4*macro(3,1)));
    alphaEigNum(kk,40)=alphaEigInf(kk,40)-...
        bndMeans(1,39)/h^2*(1-refPhi40(kk,39)/refPhi40(kk,40))-...
        2*macro(3,40)/(h*(h+4*macro(3,40)));
end

errOne=(one9(timeInds,:)-ref409)./ref409;
errTwo=(two9(timeInds,:)-ref409)./ref409;

sqrtNumPlots=4;
cellInds=round(linspace(1,39,sqrtNumPlots^2));

figure(1)
clf
buffer=.025;
plotDim=(1-buffer*(sqrtNumPlots))/sqrtNumPlots;
for i=1:sqrtNumPlots
    for j=1:sqrtNumPlots
        thisCell=cellInds((i-1)*sqrtNumPlots+j);
        axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlots-i)*(buffer+plotDim) plotDim plotDim])
        [AX H1 H2]=plotyy(time(timeInds),[errOne(:,thisCell) errTwo(:,thisCell)],time(timeInds),alphaEigInf(:,thisCell));
        set(AX(1),'YTick',[]); set(AX(1),'XTick',[]);
        set(AX(2),'YTick',[]); set(AX(2),'XTick',[]);
        set(AX(1),'XLim',[0 max(time)]); set(AX(2),'XLim',[0 max(time)]);
        text(.1,.9,num2str(thisCell),'units','normalized','FontSize',14)
    end
end

figure(2)
clf
buffer=.025;
plotDim=(1-buffer*(sqrtNumPlots))/sqrtNumPlots;
for i=1:sqrtNumPlots
    for j=1:sqrtNumPlots
        thisCell=cellInds((i-1)*sqrtNumPlots+j);
        axes('Position',[(j-1)*plotDim+buffer*(j-.5) buffer/2+(sqrtNumPlots-i)*(buffer+plotDim) plotDim plotDim])
        [AX H1 H2]=plotyy(time(timeInds),[errOne(:,thisCell) errTwo(:,thisCell)],time(timeInds),alphaEigNum(:,thisCell));
        set(AX(1),'YTick',[]); set(AX(1),'XTick',[]);
        set(AX(2),'YTick',[]); set(AX(2),'XTick',[]);
        set(AX(1),'XLim',[0 max(time)]); set(AX(2),'XLim',[0 max(time)]);
        text(.1,.9,num2str(thisCell),'units','normalized','FontSize',14)
    end
end
kill=12;