close all; 

if ex==1

    Grid.Nx=8; Grid.hx=1/Grid.Nx;
    Grid.Ny=8; Grid.hy=1/Grid.Ny;
    Grid.Nz=1; Grid.hz=1/Grid.Nz;
    Grid.K=ones(3,Grid.Nx,Grid.Ny);
    N=Grid.Nx*Grid.Ny*Grid.Nz;
    q=zeros(N,1); q([1 N])=[1 -1];
    P=TPFA(Grid,Grid.K,q);    
    contourf(P,20); axis square; colorbar; title('Pressure');

elseif ex==2
    
    Grid.Nx=32; Grid.hx=1/Grid.Nx;
    Grid.Ny=32; Grid.hy=1/Grid.Ny;
    Grid.Nz=1;  Grid.hz=1/Grid.Nz;
    Grid.K=exp(5*smooth3(smooth3(randn(3,Grid.Nx,Grid.Ny))));
    
    figure(1);    
    pcolor(log10(squeeze(Grid.K(1,:,:)))); shading flat; axis square;
    colorbar; title('log10(K)');
    
    N=Grid.Nx*Grid.Ny*Grid.Nz;
    q=zeros(N,1); q([1 N])=[1 -1];
    P=TPFA(Grid,Grid.K,q);
    
    figure(2);;
    contourf(P,20); axis square;
    colorbar; title('Pressure');
    
end
