close all;
Grid.Nx=48;  Grid.hx=20*.3048;                       % Dimension in x-direction
Grid.Ny=48; Grid.hy=10*.3048;                       % Dimension in y-direction
Grid.Nz=5;   Grid.hz=2*.3048;                        % Dimension in z-direction
Grid.N=Grid.Nx*Grid.Ny*Grid.Nz;                           % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz;                      % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;                        % Viscosities
Fluid.swc=0.2; Fluid.sor=0.2;                        % Irreducible saturations
St = 5;                                              % Maximum saturation time step
Pt = 100;                                            % Pressure time step
ND = 2000;                                           % Number of days in simulation

Q=zeros(Grid.Nx,Grid.Ny,Grid.Nz);                    % Source term for injection
IR=795*(Grid.Nx*Grid.Ny*Grid.Nz/ (60*220*85));               %   and production. Total
%Q(1,1,1)=IR; Q(Grid.Nx,Grid.Ny,Grid.Nz)=-IR; Q=Q(:); %   rate scaled to one layer
Q = Q(:);                                             % Reshape to vector

mu = 0.0;                                                   % Mean of normal
sigma = 1.0;                                                % Std of normal
Q = init_flw_trnc_norm_xin_pt_out(Grid, IR, mu, sigma, Q);  % Initialize Normal flow

load ../data/spe10;                                  % Load data from SPE 10
Grid.K=KU(:,1:Grid.Nx,1:Grid.Ny,1:Grid.Nz);          %    permeability, Layer 1
Por=pU(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz);               %    preprocessed porosity, Layer 1
Grid.por=max(Por(:),1e-3);

S=Fluid.swc*ones(Grid.N,1);                               % Initial saturation
Pc=[0; 1]; Tt=0;                                     % For production curves
for tp=1:ND/Pt;
  [P,V]=Pres(Grid,S,Fluid,Q);                        % Pressure solver
  for ts=1:Pt/St;
    S=NewtRaph(Grid,S,Fluid,V,Q,St);                 % Implicit saturation solver
    subplot('position',[0.05 .1 .4 .8]);             % Make left subplot
        pcolor(reshape(S(Grid.Nx*Grid.Ny+1:2*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');
        colorbar
        shading flat
        caxis([0 0.5]) 
                                                                  % Plot saturation
    subplot('position',[0.55 .1 .4 .8]);
    pcolor(reshape(S(2*Grid.Nx*Grid.Ny+1:3*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');     % Plot saturation
    colorbar
    caxis([0 0.5]) 
    shading flat; %caxis([Fluid.swc 1-Fluid.sor]);    %

    [Mw,Mo]=RelPerm(S(Grid.N),Fluid); Mt=Mw+Mo;           % Mobilities in well-block
    Tt=[Tt,(tp-1)*Pt+ts*St];                         % Compute simulation time
    Pc=[Pc,[Mw/Mt; Mo/Mt]];                          % Append production data
    %subplot('position',[0.55 .1 .4 .8]);             % Make right subplot
    %plot(Tt,Pc(1,:),Tt,Pc(2,:));                     % Plot production data
    %axis([0,ND,-0.05,1.05]);                         % Set correct axis
    %legend('Water cut','Oil cut');                   % Set legend
    drawnow;                                         % Force update of plot
  end
end