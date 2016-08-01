function Q = init_flw_trnc_norm_xin_pt_out(Grid, ir, mu, sigma, Q)
%initialize the total mass to 0
mass = 0.0;

q_x = zeros(1, Grid.Nx);

% Note that the portion of the  Standard Normal distribution between
% -3sigma/2 to 3sigma/2 is assumed to fit the 1..Nx
for i=1:Grid.Nx
% get the real x coordinate
  x = -1.5d0 + (i-1)*3.0d0/(Grid.Nx-1);
% Now use mu and sigma to find the pdf value at x
  result1 = sqrt(2.0d0*pi);
  arg1 = -(((x-mu)/sigma)^2.0d0/2.0d0);
  pdf = 1.0d0/(sigma*result1)*exp(arg1);
% set the value at the index equal to the pdf value at that point
  q_x(i) = pdf;
% increment the mass by the value of the pdf
  mass = mass + pdf;
end

% now rescale all the entities
q_x = q_x/mass*ir;

f1 = figure();
plot(q_x);
savefig(f1, 'inflow');

% Assign Q_x to Q
j = 1;
for i=1:Grid.Ny:Grid.Nx*Grid.Ny
  Q(i) = q_x(j);
  j = j + 1;
end

% now set the output
Q(Grid.N) = -ir;
end