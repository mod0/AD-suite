function oneD_Wave(nCell,nSec,saveEvery,filename,tsInHours)


% create subdirectory if needed to store model samples from MCMC
Files = dir(pwd); % get current directory
% check to see if results directory exists
if ~any(strcmp({Files.name},filename))
  mkdir(filename);
end

%%%%
%   Height and Radius of reactor and vector of x coordinates
%%%%
L=400;  %cm
xvec=linspace(0+L/(2*nCell),L-L/(2*nCell),nCell);
R=150;
h=L/nCell;

%Target power is the targeted fission rate in the reactor or
%sum(phi*N_9*sigma_{f,9}*h)
targetPowerDensity=100;  %W/cm^3
joulesPerFission=3.204e-11;     %3.204e-11 = J/fission
%The target fission rate given the volume of the reactor
targetFissionRate=targetPowerDensity/joulesPerFission*(L*pi*R^2);    

%Specify what portion of the domain contains the initial reaction region
nRxCell=round(.5*nCell);

%Initialize the nuclide densities and generate the matrix of macroscopic
%cross sections
%This call loads the hard-coded microscopic cross sections
getMicros;
%The nVec matrix has 3 rows and nCell columns.  
    %Row 1: Pu-239 density in each cell
    %Row 2: U-238 density in each cell
    %Row 3: FP Density in each cell
nVec=zeros(3,nCell);
nVec(1,1:nRxCell)=8e20;
nVec(1,nRxCell+1:end)=6e20;
nVec(2,1:nRxCell)=6e20;
nVec(2,nRxCell+1:end)=1.5e21;
nVec(3,:)=0;
%This call computes the macroscopic cross sections given the current nVec
    %macro - a matrix with 3 rows and nCell columns
        %Row 1: macroscopic fission cross section
        %Row 2: macroscopic absorption cross section
        %Row 3: diffusion coefficient (with H2O corrector)
    %bndMeans - a matrix with 2 rows and nCell-1 columns (interior bndrys)
        %Row 1: the harmonic diffusion cell-boundary diffusion coefficient
        %Row 2: average cell-width
getMacros;

%Calls a subroutine which gives the critical value for Sigma_ext given the
%current nuclide densities
[Sax Phi]=getSaAndFlux(macro,bndMeans);
if sum(Phi)<0
    Phi=Phi*-1;
end
%The next couple of lines fix up small negative flux values that the eigensolver
%produces
Phi(find(abs(Phi)<1e-12*max(Phi)))=abs(Phi(find(abs(Phi)<1e-12*max(Phi))));
if sum(Phi<0)>0
    error('Negative Phis!')
end
const=targetFissionRate/(sum(nVec(1,:).*h*pi*R^2*msf9.*Phi'));
PhiOld=const*Phi;

%Timestep over which phi is assumed to be linear
timeStep=60*60*tsInHours;
writeEvery=saveEvery*10;
eps=1.7e-9;
maxInner=35;
time=0;
outerCounter=0;
saveCounter=0;
warningCount=0;
%Save initial values of all variables
saveCounter=saveCounter+1;
timeStore(saveCounter)=time;
phiStore(saveCounter,:)=PhiOld;
SaxStore(saveCounter)=Sax;
nineStore(saveCounter,:)=nVec(1,:);
eightStore(saveCounter,:)=nVec(2,:);
zeroStore(saveCounter,:)=nVec(3,:);
while time<nSec
    time=time+timeStep;
    if time>=nSec
        dt=time-nSec;
        time=nSec;
        timeStep=timeStep-dt;
    end
    fprintf('\n Time (d): %.1f,',time/(60*60*24));
    outerCounter=outerCounter+1;
    nVec0=nVec;
    PhiNew=PhiOld;
    innerCounter=0;
    err=10;
    prevErr='None Yet';
    while ((err>eps) && (innerCounter<=maxInner))
        innerCounter=innerCounter+1;
        fprintf('\n        Iteration %d: Cell: ',innerCounter);
        for cl=1:nCell
            if rem(cl,round(nCell*.25))==0
%                 fprintf(['Timestep %d, Time=%.1f, Iteration %d, Updating Cell %d, Prev. Error: ' prevErr '\n'],outerCounter,time, innerCounter,cl);
                fprintf('%d...',cl);
            end
            nVec(:,cl)=computeNs(nVec0(:,cl),[PhiOld(cl) PhiNew(cl)],timeStep,1);
        end
        getMacros;
        [SaxTest Phi]=getSaAndFlux(macro,bndMeans);
        if sum(Phi)<0
            Phi=Phi*-1;
        end
        %The next couple of lines fix up small negative flux values that the eigensolver
        %produces
        Phi(find(abs(Phi)<1e-12*max(Phi)))=abs(Phi(find(abs(Phi)<1e-12*max(Phi))));
        if sum(Phi<0)>0
            error('Negative Phis!')
        end
        const=targetFissionRate/(sum(nVec(1,:).*h*pi*R^2*msf9.*Phi'));
        PhiTest=const*Phi;
%         const=targetPower/(sum(nVec(1,:).*h*msf9.*Phi'));
%         PhiTest=const*Phi;
%         PhiTest(abs(PhiTest)<1e-12*max(PhiTest))=abs(PhiTest(abs(PhiTest)<1e-12*max(PhiTest)));
%         if sum(PhiTest<0)>0
%             warning(['Negative Phis! Timestep ' num2str(outerCounter)])
%         end
%         err=max(abs((PhiNew-PhiTest)./max(PhiTest)));
        errVec=abs(PhiNew-PhiTest);
        errInds=find(PhiTest>.001*max(PhiTest));
        if length(errInds)<.25*nCell
            warning('Cells contributing to iterative convergence below 25% of problem');
            warningCount=warningCount+1;
        end
        [err errInd]=max(errVec(errInds)./PhiTest(errInds));
%         [err errInd]
        fprintf(' Error: %.2e, Warnings %d, Ext. Abs: %f',err,warningCount,SaxTest);
        prevErr=num2str(err);
        if innerCounter>maxInner
            warningCount=warningCount+1;
            warning(['Timestep ' num2str(outerCounter) ' did not Converge Phi'])
        end
        PhiNew=PhiTest;
        
    end
    PhiOld=PhiNew;
    if SaxTest<0
        fprintf('\nExternal Absorber Went Negative at t=%f years.  Ending Simulation Now.\n',...
            time/(60*60*24*365))
        saveCounter=saveCounter+1;
        timeStore(saveCounter)=time;
        phiStore(saveCounter,:)=PhiNew;
        SaxStore(saveCounter)=SaxTest;
        nineStore(saveCounter,:)=nVec(1,:);
        eightStore(saveCounter,:)=nVec(2,:);
        zeroStore(saveCounter,:)=nVec(3,:);
        break
    end
        
    if rem(outerCounter,saveEvery)==0
        saveCounter=saveCounter+1;
        timeStore(saveCounter)=time;
        phiStore(saveCounter,:)=PhiNew;
        SaxStore(saveCounter)=SaxTest;
        nineStore(saveCounter,:)=nVec(1,:);
        eightStore(saveCounter,:)=nVec(2,:);
        zeroStore(saveCounter,:)=nVec(3,:);
    end
    if rem(outerCounter,writeEvery)==0
        cd(filename)
        save phiStore phiStore;
        save nineStore nineStore;
        save eightStore eightStore;
        save zeroStore zeroStore;
        save SaxStore SaxStore;
        save timeStore timeStore;
        cd ..
    end
end
if rem(outerCounter,saveEvery)~=0
    saveCounter=saveCounter+1;
    timeStore(saveCounter)=time;
    phiStore(saveCounter,:)=PhiNew;
    SaxStore(saveCounter)=SaxTest;
    nineStore(saveCounter,:)=nVec(1,:);
    eightStore(saveCounter,:)=nVec(2,:);
    zeroStore(saveCounter,:)=nVec(3,:);
end
cd(filename)
save phiStore phiStore;
save nineStore nineStore;
save eightStore eightStore;
save zeroStore zeroStore;
save SaxStore SaxStore;
save timeStore timeStore;
cd ..
fprintf('\n')

% doPlots

% eValOld
% sxOld=Sax;
% Sax=1;
% sxNew=Sax;
% [eValNew eVec]=generalPowerIteration(1,1,1,1,1,h);
% %Uses a linear extrapolation to the next guess
% disp('Initial Criticality Search')
% while abs(eValNew-1.0)>1e-4
%     slope=(eValNew-eValOld)/(sxNew-sxOld);
%     b=eValNew-slope*sxNew;
%     Sax=(1-b)/slope;
%     sxOld=sxNew;
%     sxNew=Sax;
%     eValOld=eValNew;
%     [eValNew eVec]=generalPowerIteration(2.2*Sf,Sa,D,Dhalf,Sax,h);
% end
% eVec(1:nRxCell)=sin(pi*xvec(1:nRxCell)/xvec(nRxCell));
% eVec(nRxCell+1:nCell)=0;
% disp('Done')
% targetPower=1e10;
% pConstant=targetPower/sum(eVec.*Sf*h);
% y0(1:nCell)=eVec*pConstant;
% Sax
% 
% 
% xvec=xvec+2*D(1);
% % fluxShape=sin(xvec*pi/(400+4*D(1)));
% % const=targetPower/sum(h*fluxShape.*Sf');
% % y0(1:nCell)=const*fluxShape;
% 
% % kWOsx=2.2*Sf(1)/(Sa(1)+D(1)*(2.405*pi/(L*R))^2)
% % Sax0=2.2*Sf(1)-Sa(1)-D(1)*(2.405*pi/(L*R))^2
% 
% 
% 
% [t y]=ode15s(@yprimeWrapper,[0 nSec],y0,options);
% 
% for i=1:length(t)
%     powVal(i)=h*sum(y(i,1:nCell).*y(i,nCell+1:2*nCell)*msf9);
% end
% 
% 
% close('all')
% figure(1)
% axes('FontSize',14)
% surf(xvec,t,y(:,1:nCell))
% xlabel('X')
% ylabel('Time')
% zlabel('Flux')
% 
% figure(2)
% axes('FontSize',14)
% surf(xvec,t,y(:,nCell+1:2*nCell))
% xlabel('X')
% ylabel('Time')
% zlabel('Pu')
% 
% figure(3)
% axes('FontSize',14)
% surf(xvec,t,y(:,2*nCell+1:3*nCell))
% xlabel('X')
% ylabel('Time')
% zlabel('U')
% 
% figure(4)
% axes('FontSize',14)
% surf(xvec,t,y(:,3*nCell+1:4*nCell))
% xlabel('X')
% ylabel('Time')
% zlabel('FP')
% 
% figure(5)
% axes('FontSize',14)
% plot(xvec,y(end,1:nCell)-y(1,1:nCell))
% xlabel('X')
% ylabel('\Phi(end) - \Phi(0)');
% 
% figure(6)
% axes('FontSize',14)
% plot(xvec,y(end,nCell+1:2*nCell)-y(1,nCell+1:2*nCell))
% xlabel('X')
% ylabel('Pu(end) - Pu(0)');
% 
% figure(7)
% axes('FontSize',14)
% plot(xvec,y(end,2*nCell+1:3*nCell)-y(1,2*nCell+1:3*nCell))
% xlabel('X')
% ylabel('U(end) - U(0)');
% 
% figure(8)
% axes('FontSize',14)
% plot(xvec,y(end,3*nCell+1:4*nCell)-y(1,3*nCell+1:4*nCell))
% xlabel('X')
% ylabel('FP(end) - FP(0)');
% 
% figure(9)
% plot(t,powVal,'o')
% xlabel('Time')
% ylabel('Power')
