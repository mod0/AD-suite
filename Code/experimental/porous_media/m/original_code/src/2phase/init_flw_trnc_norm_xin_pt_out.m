function Q = init_flw_trnc_norm_xin_pt_out(Grid, ir, mu, sigma, Q)
%initialize the total mass to 0
mass = 0.0;
q_x = zeros(1, Grid.Nx);

assert(isequal(length(mu),length(sigma)));

    
for i=1:Grid.Nx
    x = (i-1.0)*2.0/(double(Grid.Nx)-1.0) - 1.0;
    pdf = 0.0;
    for j = 1:length(mu)
          arg1 = -(((x-mu(j)).^2/sigma(j))^2.0/2.0);
          pdf = pdf + 1.0/(sqrt(2.0*pi).*sigma(j))*exp(arg1);
    end
    q_x(i) = pdf;
    mass = mass + pdf;
end
% now rescale all the entities
q_x = q_x/mass*ir;

% Assign Q_x to Q
j = 1;
for i=1:Grid.Ny:Grid.Nx*Grid.Ny
  Q(i) = q_x(j);
  j = j + 1;
end

% now set the output
Q(Grid.N) = -ir;
end