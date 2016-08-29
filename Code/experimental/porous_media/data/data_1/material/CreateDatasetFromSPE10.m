function CreateDatasetFromSPE10(filename, nx, ny, nz)
    data = load(filename, 'pU', 'KU');             
    
    % expects variables pU and KU upon loading file
    maxnx = size(data.pU, 1);
    maxny = size(data.pU, 2);
    maxnz = size(data.pU, 3);
    
    % ratio of maxnx to nx, maxny to ny, maxnz to nz
    rnx = int32(floor(maxnx/nx));
    rny = int32(floor(maxny/ny));
    rnz = int32(floor(maxnz/nz));
    
    if(nx > maxnx || ny > maxny || nz > maxnz)
        error('Input dimensions exceed dimension of dataset'); 
    end

    % Compute hX, hY, hZ
    hX = double(1200/nx) * 0.3048;
    hY = double(2200/ny) * 0.3048;
    hZ = double(170/nz)  * 0.3048;

    % Compute volume
    vol = hX * hY * hZ;

    % Compute ir_const
    ir_const = 795;
    
    pU_new = zeros(nx, ny, nz);
    KU_new = zeros(3, nx, ny, nz);
    
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                pU_temp = data.pU((i-1) * rnx + 1:i*rnx, ...
                                  (j-1) * rny + 1:j*rny, ...
                                  (k-1) * rnz + 1:k*rnz);
                pU_new(i, j, k) = mean(pU_temp(:));
                
                KU_temp = data.KU(1, (i-1) * rnx + 1:i*rnx, ...
                                  (j-1) * rny + 1:j*rny, ...
                                  (k-1) * rnz + 1:k*rnz);
                KU_new(1,i,j,k) = mean(KU_temp(:));
                
                KU_temp = data.KU(2, (i-1) * rnx + 1:i*rnx, ...
                                  (j-1) * rny + 1:j*rny, ...
                                  (k-1) * rnz + 1:k*rnz);
                KU_new(2,i,j,k) = mean(KU_temp(:));
                
                KU_temp = data.KU(3, (i-1) * rnx + 1:i*rnx, ...
                                  (j-1) * rny + 1:j*rny, ...
                                  (k-1) * rnz + 1:k*rnz);
                KU_new(3,i,j,k) = mean(KU_temp(:));
            end
        end
    end
    
    % write porosity data
    ncid = netcdf.create('porosity.nc','CLOBBER');
    ndim = netcdf.defDim(ncid, 'Ndim', nx*ny*nz);
    pUvarid = netcdf.defVar(ncid, 'porosity', 'double', ndim);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid, pUvarid, pU_new(:));
    netcdf.close(ncid);
    
    % write permeability data
    ncid = netcdf.create('permeability.nc','CLOBBER');
    fixdim = netcdf.defDim(ncid, 'FixedDim', 3);
    nxdim = netcdf.defDim(ncid, 'Nx', nx);
    nydim = netcdf.defDim(ncid, 'Ny', ny);
    nzdim = netcdf.defDim(ncid, 'Nz', nz);
    KUvarid = netcdf.defVar(ncid, 'permeability', 'double', [fixdim, nxdim, nydim, nzdim]);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid, KUvarid, KU_new);
    netcdf.close(ncid);
end