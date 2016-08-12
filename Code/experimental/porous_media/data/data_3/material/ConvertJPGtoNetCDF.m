function [ ] = ConvertJPGtoNetCDF( file1, file2, file3 )
%This function reads an jpg file and returns the hsv img

if nargin < 1 
    file1 = 'porosity_large.jpg';
    file2 = 'permeability_large.jpg';
    file3 = 'saturation_large.jpg';
end

tmp1 = rgb2hsv(imread(file1));
tmp2 = rgb2hsv(imread(file2));
tmp3 = rgb2hsv(imread(file3));
tmp1 = 1.0-tmp1;
tmp2 = 1.0-tmp2;
tmp3 = 1.0-tmp3;

assert(isequal(size(tmp1),size(tmp2)))
assert(isequal(size(tmp2),size(tmp3)))

% This needs to be corrected! 
por(:,:,1) = tmp1(:,:,3);
por(:,:,2) = tmp1(:,:,3);

perm(1,:,:,:) = tmp2(:,:,3);
perm(2,:,:,:) = tmp2(:,:,3);
perm(3,:,:,:) = tmp2(:,:,3);
perm(:,:,:,1) = perm;
perm(:,:,:,2) = perm;

sat(:,:,1) = tmp3(:,:,3);
sat(:,:,2) = tmp3(:,:,3);


% Set some parameters

Nf = 3;
NX = size(tmp1(:,:,1),1);
NY = size(tmp1(:,:,1),2);
NZ = 2;
hX = 1.0/NX;
hY = 1.0/NY;
hZ = 1.0/NZ;
vw = 1.e-04;
vo = 1.e-04;
swc = 0.0;
sor = 0.0;
V   = hX*hY*hZ;
St = 5;
Pt = 100;
ND = 2000;
ir_const = 1.0;
solver_inner = 64;
solver_outer = 100000;



%Create Parameter files
%Dimensions

ncid = netcdf.create('../parameters1.nc','NC_WRITE');
dimid_nx = netcdf.defDim(ncid,'NX',1);
dimid_ny = netcdf.defDim(ncid,'NY',1);
dimid_nz = netcdf.defDim(ncid,'NZ',1);
dimid_hx = netcdf.defDim(ncid,'hX',1);
dimid_hy = netcdf.defDim(ncid,'hY',1);
dimid_hz = netcdf.defDim(ncid,'hZ',1);
dimid_vw = netcdf.defDim(ncid,'vw',1);
dimid_vo = netcdf.defDim(ncid,'vo',1);
dimid_swc = netcdf.defDim(ncid,'swc',1);
dimid_sor = netcdf.defDim(ncid,'sor',1);
dimid_V = netcdf.defDim(ncid,'V',1);
dimid_st = netcdf.defDim(ncid,'St',1);
dimid_pt = netcdf.defDim(ncid,'Pt',1);
dimid_nd = netcdf.defDim(ncid,'ND',1);
dimid_ir = netcdf.defDim(ncid,'ir_const',1);
dimid_inner = netcdf.defDim(ncid,'solver_inner',1);
dimid_outer = netcdf.defDim(ncid,'solver_outer',1);

varid_nx = netcdf.defVar(ncid,'NX','int',dimid_nx);
varid_ny = netcdf.defVar(ncid,'NY','int',dimid_ny);
varid_nz = netcdf.defVar(ncid,'NZ','int',dimid_nz);
varid_hx = netcdf.defVar(ncid,'hX','double',dimid_hx);
varid_hy = netcdf.defVar(ncid,'hY','double',dimid_hy);
varid_hz = netcdf.defVar(ncid,'hZ','double',dimid_hz);
varid_vw = netcdf.defVar(ncid,'vw','double',dimid_vw);
varid_vo = netcdf.defVar(ncid,'vo','double',dimid_vo);
varid_swc = netcdf.defVar(ncid,'swc','double',dimid_swc);
varid_sor = netcdf.defVar(ncid,'sor','double',dimid_sor);
varid_V = netcdf.defVar(ncid,'V','double',dimid_V);
varid_st = netcdf.defVar(ncid,'St','int',dimid_st);
varid_pt = netcdf.defVar(ncid,'Pt','int',dimid_pt);
varid_nd = netcdf.defVar(ncid,'ND','int',dimid_nd);
varid_ir = netcdf.defVar(ncid,'ir_const','double',dimid_ir);
varid_inner = netcdf.defVar(ncid,'solver_inner','int',dimid_inner);
varid_outer = netcdf.defVar(ncid,'solver_outer','int',dimid_outer);

netcdf.endDef(ncid)
netcdf.putVar(ncid,varid_nx,NX);
netcdf.putVar(ncid,varid_ny,NY);
netcdf.putVar(ncid,varid_nz,NZ);
netcdf.putVar(ncid,varid_hx,hX);
netcdf.putVar(ncid,varid_hy,hY);
netcdf.putVar(ncid,varid_hz,hZ);
netcdf.putVar(ncid,varid_vw,vw);
netcdf.putVar(ncid,varid_vo,vo);
netcdf.putVar(ncid,varid_swc,swc);
netcdf.putVar(ncid,varid_sor,sor);
netcdf.putVar(ncid,varid_V,V);
netcdf.putVar(ncid,varid_st,St);
netcdf.putVar(ncid,varid_pt,Pt);
netcdf.putVar(ncid,varid_nd,ND);
netcdf.putVar(ncid,varid_ir,ir_const);
netcdf.putVar(ncid,varid_inner,solver_inner);
netcdf.putVar(ncid,varid_outer,solver_outer);
netcdf.close(ncid) 


%Permeability
ncid = netcdf.create('../parameters2.nc','NC_WRITE');
dimidfix = netcdf.defDim(ncid,'FixedDim',Nf);
dimidnx  = netcdf.defDim(ncid,'NX',NX);
dimidny  = netcdf.defDim(ncid,'NY',NY);
dimidnz  = netcdf.defDim(ncid,'NZ',NZ);
varid    = netcdf.defVar(ncid,'permeability','NC_DOUBLE',[dimidfix dimidnx dimidny dimidnz]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid, perm);
netcdf.close(ncid);

%Porosity
ncid = netcdf.create('../parameters3.nc','NC_WRITE');
dimidnx = netcdf.defDim(ncid,'NX',NX);
dimidny = netcdf.defDim(ncid,'NY',NY);
dimidnz = netcdf.defDim(ncid,'NZ',NZ);
varid = netcdf.defVar(ncid,'porosity','NC_DOUBLE',[dimidnx dimidny dimidnz]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid, por);
netcdf.close(ncid);

%Saturation
ncid = netcdf.create('../parameters4.nc','NC_WRITE');
dimidnx = netcdf.defDim(ncid,'NX',NX);
dimidny = netcdf.defDim(ncid,'NY',NY);
dimidnz = netcdf.defDim(ncid,'NZ',NZ);
varid = netcdf.defVar(ncid,'saturation','NC_DOUBLE',[dimidnx dimidny dimidnz]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid, sat);
netcdf.close(ncid);



% Create file for independent variables
n_dim = 6;
x     = zeros(n_dim ,1);
x(1)  = -1.0d+0;
x(2)  =  0.0d+0;
x(3)  = +1.0d+0;
x(4)  = +1.0d-2;
x(5)  = +1.0d-2;
x(6)  = +1.0d-2;

ncid = netcdf.create('../x.nc','NC_WRITE');
dimid_ndim = netcdf.defDim(ncid,'n_dim',1);
dimid_x    = netcdf.defDim(ncid,'x',n_dim); 
varid_ndim = netcdf.defVar(ncid,'n_dim','int',dimid_ndim);
varid_x = netcdf.defVar(ncid,'x','double',dimid_x);
netcdf.endDef(ncid)
netcdf.putVar(ncid,varid_ndim,n_dim);
netcdf.putVar(ncid,varid_x,x);
netcdf.close(ncid) 


% Create file for dependent variables
m_dim = 1;
y     = zeros(m_dim ,1);

ncid = netcdf.create('../y.nc','NC_WRITE');
dimid_mdim = netcdf.defDim(ncid,'m_dim',1);
dimid_y    = netcdf.defDim(ncid,'y',m_dim); 
varid_mdim = netcdf.defVar(ncid,'m_dim','int',dimid_mdim);
varid_y = netcdf.defVar(ncid,'y','double',dimid_y);
netcdf.endDef(ncid)
netcdf.putVar(ncid,varid_mdim,m_dim);
netcdf.putVar(ncid,varid_y,y);
netcdf.close(ncid) 


% Some output
ncdisp('../parameters1.nc')
ncdisp('../parameters2.nc')
ncdisp('../parameters3.nc')
ncdisp('../parameters4.nc')
ncdisp('../x.nc')
ncdisp('../y.nc')
subplot(1,3,1)
    pcolor(tmp1(:,:,3));
    shading interp
subplot(1,3,2)
    pcolor(tmp2(:,:,3))
    shading interp
subplot(1,3,3)
    pcolor(tmp3(:,:,3))
    shading interp

end

