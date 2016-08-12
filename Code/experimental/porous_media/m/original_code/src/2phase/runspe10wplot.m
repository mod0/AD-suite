close all;
addpath('../TPFA')

Grid.Nx=48;  Grid.hx=20*.3048;                        % Dimension in x-direction
Grid.Ny=48; Grid.hy=10*.3048;                         % Dimension in y-direction
Grid.Nz=5;   Grid.hz=2*.3048;                         % Dimension in z-direction
Grid.N=Grid.Nx*Grid.Ny*Grid.Nz;                       % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz;                       % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;                         % Viscosities
Fluid.swc=0.2; Fluid.sor=0.2;                         % Irreducible saturations
St = 5;                                               % Maximum saturation time step
Pt = 100;                                             % Pressure time step
ND = 2000;                                            % Number of days in simulation

IR=795*(Grid.Nx*Grid.Ny*Grid.Nz/ (60*220*85));        % and production. Total

mu    = [ -1 0 1];                                    % Mean of normal
sigma = 0.01*ones(length(mu));                        % Std of normal

Q = zeros(Grid.Nx,Grid.Ny,Grid.Nz);                     % Source term for injection
Q = init_flw_trnc_norm_xin_pt_out(Grid, IR, mu, sigma, Q);  % Initialize Normal flow
Q = Q(:);                                             % Reshape to vector

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
    [Mw,Mo]=RelPerm(S(Grid.N),Fluid);                % Mobilities in well-block
    Mt=Mw+Mo; 
    Tt=[Tt,(tp-1)*Pt+ts*St];                         % Compute simulation time
    Pc=[Pc,[Mw/Mt; Mo/Mt]];                          % Append production data
    
    subplot(2,3,1)
        pcolor(reshape(S(Grid.Nx*Grid.Ny+1:2*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');   % Plot saturation
        colorbar
        caxis([0 1]) 
        shading flat
                                                   
    subplot(2,3,2)
        pcolor(reshape(S(2*Grid.Nx*Grid.Ny+1:3*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');     % Plot saturation
        colorbar
        caxis([0 1]) 
        shading flat

    subplot(2,3,3)
        pcolor(reshape(S(3*Grid.Nx*Grid.Ny+1:4*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');     % Plot saturation
        colorbar
        caxis([0 1]) 
        shading flat

    subplot(2,3,4)
        pcolor(reshape(S(4*Grid.Nx*Grid.Ny+1:5*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');     % Plot saturation
        colorbar
        caxis([0 1]) 
        shading flat

    subplot(2,3,4)
        pcolor(reshape(S(4*Grid.Nx*Grid.Ny+1:5*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');     % Plot saturation
        colorbar
        caxis([0 1]) 
        shading flat

    subplot(2,3,5)
        plot(Tt,Pc(1,:),Tt,Pc(2,:));                     % Plot production data
        axis([0,ND,-0.05,1.05]);                         % Set correct axis
        legend('Water cut','Oil cut');                   % Set legend
    
    drawnow;                                             % Force update of plot
  end
end