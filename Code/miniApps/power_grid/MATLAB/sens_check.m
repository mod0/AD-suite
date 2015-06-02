% Script to check sensitivities
function [costfunction,adjsens] = sens_check(pm,tend)
global  pmax omegaB omegaS H D thetaS t0 tf tcl beta c
global g dg dt adjsens tempder epsilon


pmax = 2.0877;
omegaB = 120*pi;
omegaS = 1;
H = 5.0;
D = 5.0;
thetaS = 1.0;
beta = 2;
c=10000;
t0 = 0;
tend = 10.0;
tf = 0.1;   % fault-on time
tcl = 0.2; % fault-off time
dt = 0.01;

perturb = 1.0; %% initial condition perturbation (this has nothing to do with any numerical scheme)

epsilon = 1e-8; %% Finite difference perturbation

g = 0;
dg = 0;
adjsens = 0;
tempder = 0;


int_method = @CrankNicholson;
%pm = 1.2;
pm
xinit = [asin(perturb*pm/pmax);perturb*1.0]; %% Perturb initial conditions so we can see a transient
tspan = t0:dt:tend;

mass = eye(2);
options = odeset('Mass',mass,'MassSingular','no','Jacobian',@(t,x)ForwardOdeRHSJacobian(t,x,pm),'OutputFcn',@ForwardOdeOutputFcn);

[tout, yout] = int_method(@(t,x)ForwardOdeRhs(t,x,pm),tspan,xinit,options);

figure(1);
plot(tout,yout(:,1))
hold on
plot(tout,thetaS);

% figure(2);
% plot(yout(:,1),yout(:,2));
% xlabel('delta (rad)');
% ylabel('omega (pu)');

%% Sensitivity calculation using finite differencing
%% Uncomment this part when you want to check the accuracy of the gradient calculated
%% by adjoints. dg_dp is the gradient
%g1 = g;
%g = 0;
%dg = 0;
%pm1 = pm + epsilon;
%xinit = [asin(pm1/pmax);1.0];
%tspan = t0:dt:tend;
%[tout,yout] = int_method(@(t,x)ForwardOdeRhs(t,x,pm1),tspan,xinit,options);
%dg_dp = -1 + (g-g1)/epsilon;

%% Sensitivitity calculation using adjoints
xinit = [0;0];
tspan = fliplr(tspan);
options = odeset('Mass',mass,'MassSingular','no','Jacobian',@(t,x)AdjointOdeRHSJacobian(t,x,pm,tout,yout),'OutputFcn',@(t,x,flag)AdjointOdeOutputFcn(t,x,flag,pm,tout,yout,perturb));
[tb, yb] = int_method(@(t,x)AdjointOdeRhs(t,x,pm,tout,yout),tspan,xinit,options);
%%%%%%%


% figure(3);
% plot(tb,yb);
%
% figure(4);
% plot(yb(:,1),yb(:,2));
% xlabel('lambda_1');
% ylabel('lambda_2');

adjsens
costfunction=-pm+g;

end



function [tout,yout] = CrankNicholson(odefun,tspan,x0,options)

tout = zeros(length(tspan),1);
yout = zeros(length(tspan),length(x0));
tout(1) = tspan(1);
yout(1,:) = x0';

h = tspan(2) - tspan(1);

t = tout(1);

terminate = feval(options.OutputFcn,t,x0,'init');

max_it = 100;


df0 = feval(odefun,t,x0);
for i=2:length(tspan)
    x = x0;
    t = tspan(i);
    
    converged = 0;
    
    df = feval(odefun,t,x);
    f = x - x0 - 0.5*h*(df+df0);
    if (norm(f,inf) < 1e-8)
        converged = 1;
        x0 = x;
        df0 = df;
    end
    ctr = 0;
    while(~converged && ctr < max_it)
        ctr = ctr + 1;
        
        J = feval(options.Jacobian,t,x);
        
        J = eye(length(x)) - 0.5*h*J;
        
        dx = -(J\f);
        x  = x + dx;
        
        df = feval(odefun,t,x);
        f = x - x0 - 0.5*h*(df+df0);
        if (norm(f,inf) < 1e-8)
            converged = 1;
            x0 = x;
            df0 = df;
        end
    end
    
    terminate = feval(options.OutputFcn,t,x,'');
    
    
    tout(i) = t;
    yout(i,:) = x';
end

terminate = feval(options.OutputFcn,t,x,'done');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward Integration routines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ForwardOdeRhs(t,x,pm)
global omegaB omegaS pmax D H tf tcl

y=zeros(2,1);
delta = x(1);
omega = x(2);

if(t > tf && t <= tcl)
    pmax1 = 0;
else
    pmax1 = pmax;
end
y(1) = omegaB*(omega - omegaS);
y(2) = omegaS*(pm - pmax1*sin(delta) -D*(omega - omegaS))/(2*H);

end

function [J] = ForwardOdeRHSJacobian(t,x,pm)
global omegaB omegaS pmax D H tf tcl
y=zeros(2,1);
delta = x(1);
omega = x(2);

if(t > tf && t <= tcl)
    pmax1 = 0;
else
    pmax1 = pmax;
end

J = [0                                omegaB;
    -omegaS*pmax1*cos(delta)/(2*H)    -D/(2*H)];

end

function [terminate] = ForwardOdeOutputFcn(t,x,flag)
terminate = 0;

global c thetaS g dg dt beta

if isempty(flag)
    
    delta = x(1);
    
    g = g + 0.5*dt*(dg+c*max(0,(delta-thetaS))^beta);
    
    dg = c*max(0,(delta-thetaS))^beta;
    
elseif flag == 'init'
    delta = x(1);
    g = 0;
    dg = c*max(0,(delta-thetaS))^beta;
    return;
elseif flag == 'done'
    g
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adjoint Integration routines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = AdjointOdeRhs(t,x,pm,tout,yout)
global thetaS beta c
y=zeros(2,1);

idx = find(t == tout);

J = ForwardOdeRHSJacobian(tout(idx),yout(idx,:)',pm);

delta = yout(idx,1);
forcing = [c*beta*max(0,delta - thetaS)^(beta-1);0];

y = -J'*x - forcing;

end

function [J] = AdjointOdeRHSJacobian(t,x,pm,tout,yout)
idx = find(t == tout);

J = ForwardOdeRHSJacobian(tout(idx),yout(idx,:)',pm);

J = -J';

end


function [terminate] = AdjointOdeOutputFcn(t,x,flag,pm,tout,yout,perturb)
terminate = 0;

global thetaS g dg dt H adjsens tempder pmax epsilon beta c

idx = t == tout;
if isempty(flag)
    
    adjsens = adjsens + 0.5*dt*(tempder + [0,1/(2*H)]*x);
    
    tempder = [0,1/(2*H)]*x;
    
elseif flag == 'init'
    adjsens = 0;
    delta = yout(idx,1);
    tempder = [c*beta*max(0,delta - thetaS)^(beta-1)];
    %tempder = 0;
    return;
elseif flag == 'done'
    idx = t == tout;
    % pm = pm + epsilon;
    xinit = [asin(perturb*pm/pmax);perturb*1.0]; %% Perturb initial conditions so we can see a transient
    % dx0_dp = (xinit-yout(idx,:)')/epsilon;
    dx0_dp = [(1/pmax)*(1-(pm/pmax)^2)^(-0.5);0];
    adjsens = -1+adjsens + dx0_dp'*x;
end
end

