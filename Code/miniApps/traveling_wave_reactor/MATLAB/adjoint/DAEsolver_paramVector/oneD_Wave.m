function [t y p rxParams]=oneD_Wave(nCell,nSec,saveEveryInDays,p)
%%%%
%   ODE OPTIONS
%%%%
massMat=eye(4*nCell+1);
for i=3*nCell+1:4*nCell+1
    massMat(i,i)=0;
end
options=odeset('Mass',massMat,...
    'MStateDependence','none',...
    'MassSingular','yes',...
    'AbsTol',1e4,...
    'RelTol',1e-6,...
    'Events',@stopODE);
%     'Jacobian',@yprime_forJacobian);

%%%%
%   Height and Radius of reactor and vector of x coordinates
%%%%
rxParams.nCell=nCell;
rxParams.L=400;
rxParams.R=150;
rxParams.B2geo=(2.405/rxParams.R)^2;
rxParams.h=rxParams.L/nCell;
rxParams.nAlg=nCell+1;
rxParams.nDiff=3*nCell;
%Target power is the targeted fission rate in the reactor or
%sum(phi*N_9*sigma_{f,9}*h)
rxParams.targetPowerDensity=100;  %W/cm^3
rxParams.joulesPerFission=3.204e-11;     %3.204e-11 = J/fission
%The target fission rate given the volume of the reactor
rxParams.targetFissionRate=...
    rxParams.targetPowerDensity/rxParams.joulesPerFission*(rxParams.L*pi*rxParams.R^2);    
%Specify what portion of the domain contains the initial reaction region
rxParams.nRxCell=round(.5*nCell);

%Initialize the parameters vector
% p=[msf9 msa9 mst9 msa8 mst8 msa0 mst0 Gamma nu alphaD];
% p=[1000.4e-24 (1026-7.968)*1e-24 1026e-24 ...
%             500.717e-24 600.09e-24 ...
%             5e-24 10e-24 ...
%             .1*(1-1000.4/(1026-7.968)) 2.2 500]';
%Initialize the vector of differential uknowns
    %First nCell unknowns: N_9
    %Second nCell unknowns: N_8
    %Third nCell unknowns: N_0
    %Fourth nCell unknowns phi
    %last unknown Sigma_a^{ext}
y0=zeros(4*rxParams.nCell+1,1);
y0(1:rxParams.nRxCell)=8e20;
y0(rxParams.nRxCell+1:rxParams.nCell)=6e20;
y0(rxParams.nCell+1:rxParams.nCell+rxParams.nRxCell)=6e20;
y0(rxParams.nCell+rxParams.nRxCell+1:2*rxParams.nCell)=1.5e21;

%Calls a subroutine which gives the critical value for Sigma_ext given the
%current nuclide densities
[Sax Phi]=getSaAndFlux(y0,p,rxParams);
if sum(Phi)<0
    Phi=Phi*-1;
end
%The next couple of lines fix up small negative flux values that the eigensolver
%produces
Phi(find(abs(Phi)<1e-12*max(Phi)))=abs(Phi(find(abs(Phi)<1e-12*max(Phi))));
if sum(Phi<0)>0
    error('Negative Phis!')
end
const=rxParams.targetFissionRate/(sum(y0(1:rxParams.nCell).*rxParams.h*pi*rxParams.R^2*p(1).*Phi));
PhiOld=const*Phi;

%Now we will set up for the Matlab solver to compute the solution
%y0 is the initial condition for the system
y0(3*nCell+1:4*nCell)=PhiOld(1:end);
y0(4*nCell+1)=Sax*10^10;

%Specify the times at which we want the solution
saveEvery=saveEveryInDays*24*60*60;
tVec=[0:saveEvery:nSec];
if tVec(end)~=nSec
    tVec=[tVec nSec];
end
%OK - Call the solver
tracker=0;
save tracker tracker;
fprintf('Simulation Time (yr):');
[t y]=ode23t(@(t,y) yprime(t,y,p,rxParams),tVec,y0,options);
tyr=t/(60*60*24*365);
fprintf('\n');

frdMat=rxParams.h*pi*rxParams.R^2*p(1)*y(:,3*nCell+1:4*nCell).*y(:,1:nCell);
reactorPower=rxParams.joulesPerFission*sum(frdMat,2)*1e-6;
xvec=linspace(0+rxParams.L/(2*nCell),rxParams.L-rxParams.L/(2*nCell),nCell);

% figure(1)
% clf
% axes('FontSize',14);
% plot(tyr,reactorPower)
% xlabel('Time (yr)');
% ylabel('Reactor Power (MegaWatts)');
% 
% figure(2)
% clf
% axes('FontSize',14);
% surf(xvec,tyr,y(:,3*nCell+1:4*nCell))
% % set(gca,'zscale','log');
% xlabel('Space (cm)')
% ylabel('Time (yr)');
% zlabel('Neutron Flux (cm^{-2}s^{-1})');
% 
% figure(3)
% clf
% axes('FontSize',14);
% plot(tyr,y(:,end)/10^10);
% xlabel('Time (yr)');
% ylabel('External Absorber XS (cm^{-1})')
% 
% figure(4)
% clf
% axes('FontSize',14);
% surf(xvec,tyr,y(:,1:nCell))
% % set(gca,'zscale','log');
% xlabel('Space (cm)');
% ylabel('Time (yr)');
% zlabel('Pu-239 Density')


