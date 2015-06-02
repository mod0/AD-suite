function [sol solBack rat]=testAnalyticDAE(doPlot,tol)
%Add the dir containig my hacked ode23t solver
addpath('builtIns/')

%Mass Matrix for forward solve
M=eye(2);
M(2,2)=0;
options=odeset('Mass',M,...
    'MStateDependence','none',...
    'MassSingular','yes','relTol',tol);

%Initial conditions/parameters
y0=[1;0];
t_end=10;
p=pi;

%Call the forward solver
[sol]=ode23t_hayes(@(t,y) yprimeTest(t,y,p),[0 t_end],y0,options);

%Evaluate the solution
t=linspace(0,t_end,100);
yNum=deval(sol,t);

%Analytic Solution 
yAnl=y0(1)*exp(-p*t.^2/400);

%Analytic Derivative w.r.t. p
dydp_anl=-y0(1)*t.^2/400.*exp(-p*t.^2/400);

%Finite Difference approximation to dy/dp for verification
dydp_FD=(yAnl-y0(1)*exp(-1*(p+1e-3*p)*t.^2/400))/(-1e-3*p);

if doPlot
    %Plot function
    figure(1)
    axes('FontSize',14)
    plot(t,yAnl,t,yNum(1,:),'go')
    xlabel('Time')
    ylabel('Response')
    legend('Analytic','Numerical')
    
    %Plot dydp
    figure(2)
    axes('FontSize',14);
    plot(t,dydp_anl,t,dydp_FD,'go')
    xlabel('Time');
    ylabel('dydp')
    legend('Analytic','Numerical')
end

%backwards integration (g=Y)
%lambda0:
    %lambda0(1) -- differential adjoint unknown
    %lambda0(2) -- algebraic adjoint unknown
    %lambda0(3) -- integration of truncation term
    %lambda0(4) -- integration for dg/dp
lambda0=zeros(4,1);
lambda0(1)=1;
%algebraic constraint determines \lambda^a
lambda0(2)=-p*yNum(1,end)*lambda0(1);
lambda0(3:4)=0;

%Mass matrix for backward integration
M=eye(4);
M(2,2)=0;
options=odeset('Mass',M,...
    'MStateDependence','none',...
    'MassSingular','yes','relTol',tol/10);
[solBack]=ode23t(@(t,y) lambdaPrimeTest(t,y,sol,p),[t_end 0],lambda0,options);


lam0=deval(solBack,0);
globalError=yAnl(end)-yNum(1,end);
estGlobalErr=lam0(3);
byHandErr=estimateGlobalError(sol,solBack,25);
rat=globalError/byHandErr;

%Output
fprintf('True Global Error: %.4e\n',globalError);
fprintf('Est. Global Error Tight Integration %.4e\n',byHandErr);
fprintf('Est. Global Error: %.4e (Rel Error = %.3e)\n',estGlobalErr,(estGlobalErr-globalError)/globalError);
fprintf('True dgdp: %.4e\n',dydp_anl(end));
fprintf('Est. dpgp: %.4e (Rel Error = %.3e)\n',lam0(4),(lam0(4)-dydp_anl(end))/dydp_anl(end));
kill=12;