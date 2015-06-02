function [t y]=oneD_Wave(nCell,nSec)
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

%Now we will set up for the Matlab solver to compute the solution
%y0 is the initial condition for the system
y0=zeros(4*nCell+1,1);
y0(1:nCell)=nVec(1,:);
y0(nCell+1:2*nCell)=nVec(2,:);
y0(2*nCell+1:3*nCell)=nVec(3,:);
y0(3*nCell+1:4*nCell)=PhiOld(1:end);
y0(4*nCell+1)=Sax*10^10;

%Specify the times at which we want the solution
saveEvery=60*60*24*20;
tVec=[0:saveEvery:nSec];
if tVec(end)~=nSec
    tVec=[tVec nSec];
end
%OK - Call the solver
tracker=0;
save tracker tracker;
fprintf('Simulation Time (yr):');
[t y]=ode23t(@yprime,tVec,y0,options);
tyr=t/(60*60*24*365);
fprintf('\n\n');

frdMat=h*pi*R^2*msf9*y(:,3*nCell+1:4*nCell).*y(:,1:nCell);
reactorPower=joulesPerFission*sum(frdMat,2)*1e-6;

figure(1)
clf
axes('FontSize',14);
plot(tyr,reactorPower)
xlabel('Time (yr)');
ylabel('Reactor Power (MegaWatts)');

figure(2)
clf
axes('FontSize',14);
surf(xvec,tyr,y(:,3*nCell+1:4*nCell))
% set(gca,'zscale','log');
xlabel('Space (cm)')
ylabel('Time (yr)');
zlabel('Neutron Flux (cm^{-2}s^{-1})');

figure(3)
clf
axes('FontSize',14);
plot(tyr,y(:,end)/10^10);
xlabel('Time (yr)');
ylabel('External Absorber XS (cm^{-1})')

figure(4)
clf
axes('FontSize',14);
surf(xvec,tyr,y(:,1:nCell))
% set(gca,'zscale','log');
xlabel('Space (cm)');
ylabel('Time (yr)');
zlabel('Pu-239 Density')


