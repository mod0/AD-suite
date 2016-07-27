Grid.Nx=60;  Grid.hx=20*.3048;              % Dimension in x-direction
Grid.Ny=220; Grid.hy=10*.3048;              % Dimension in y-direction
Grid.Nz=1;   Grid.hz= 2*.3048;              % Dimension in z-direction
Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz;         % Local variables
N=Nx*Ny*Nz;                                 % Total number of grid points

Ex=(Nx-1)*Ny*Nz;                            % Number of edges in x-direction
Ey=Nx*(Ny-1)*Nz;                            % Number of edges in y-direction
Ez=Nx*Ny*(Nz-1);                            % Number of edges in z-direction
E=Ex+Ey+Ez;                                 % Total number of edges in grid

q=zeros(E+N,1);                             % Right-hand side
q(E+1)=1;                                   % Injection in block (1,1,1)
q(E+N)=-1;                                  % Production in block (Nx,Ny,Nz)

load ../data/spe10;                         % Read permeability from SPE10
Layer = 1;                                  %   which layer to extract
Grid.K=KU(:,1:Nx,1:Ny,Layer);               %   extract horizontal layer

B=GenB(Grid,Grid.K);                        % Compute B-block of matrix
C=GenC(Grid);                               % Compute C-block of matrix
A=[B,C';-C,sparse(N,N)];                    % Assemble matrix
A(E+1,E+1)=A(E+1,E+1)+1;

warning off;
x=A\q;                                      % Solve linear system
warning on;

v=x(1:E);                                   % Extract velocities
vx=reshape(v(1:Ex),Nx-1,Ny,Nz);             %   x-component
vy=reshape(v(Ex+1:E-Ez),Nx,Ny-1,Nz);        %   y-component
vz=reshape(v(E-Ez+1:E),Nx,Ny,Nz-1);         %   z-component
p=reshape(x(E+1:E+N),Nx,Ny,Nz);             % Extract pressure
