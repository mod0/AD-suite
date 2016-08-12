function Q = InitInOut(Grid, ir, mu, sigma)
%initialize the total mass to 0
Q   = zeros(Grid.N, 1);
q_x = zeros(Grid.Nx,1);
    
mass = 0.0;
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
Q(Grid.N)          = -ir/2;
Q(Grid.Ny)         = -ir/2;

end