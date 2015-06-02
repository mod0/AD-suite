%Microscopic Cross-section data.  Computed for nuclides at 0.0253eV using
%atom.kaeri.re.kr

%Uranium-238
mst8=600.09e-24;
msa8=500.717e-24;


%Plutonium-239
mst9=1026e-24;
msf9=1000.4e-24;
msa9=(1026-7.968)*1e-24;

%Representative Fission Product
mst0=10e-24;
msa0=5e-24;

%Gamma = probability that non-fission in 239-Pu will lead to a future
%absorption in a fissile nuclide
Gamma=.1*(1-msf9/msa9);

%neutron speed
nSpeed=1.5e6;

% %Uranium-238
% mst8=85e-24;
% msa8=80e-24;
% 
% 
% %Plutonium-239
% mst9=820e-24;
% msf9=700e-24;
% msa9=800e-24;
% 
% %Representative Fission Product
% mst0=1000e-24;
% msa0=500e-24;
% 
% %Gamma = probability that non-fission in 239-Pu will lead to a future
% %absorption in a fissile nuclide
% Gamma=.5*(1-msf9/msa9);
% 
% %neutron speed
% nSpeed=1.5e7;

% %Uranium-238
% mst8=12.09e-24;
% msa8=(2.717e-24+277e-24)/2;
% 
% 
% %Plutonium-239
% mst9=1026e-24;
% msf9=747.4e-24;
% msa9=(1026-7.968)*1e-24;
% 
% %Representative Fission Product
% mst0=1000e-24;
% msa0=500e-24;
% 
% %Gamma = probability that non-fission in 239-Pu will lead to a future
% %absorption in a fissile nuclide
% Gamma=.5*(1-msf9/msa9);
% 
% %neutron speed
% nSpeed=1.5e7;