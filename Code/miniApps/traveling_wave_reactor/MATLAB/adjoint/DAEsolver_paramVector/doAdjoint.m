function lambda = doAdjoint(t,y,p,rxParams,metFun)
nAlg=rxParams.nAlg;
nDiff=rxParams.nDiff;
diffRng=1:nDiff;
algRng=nDiff+1:(nAlg+nDiff);
lambda=zeros(size(y));

[g g_p g_xd g_xa]=metFun(t,y,p,rxParams);
%for fertile inventory at t_f
% g_xd=zeros(nAlg+nDiff);
% g_xd(rxParams.nCell+1:2*rxParams.nCell)=rxParams.h*pi*rxParams.R^2;

dfdx=compute_dFdx(y(end,:),p,rxParams);
ldt=g_xd'-dfdx(algRng,diffRng)'*(dfdx(algRng,algRng)'\g_xa');
lat=(dfdx(algRng,algRng)')\(-dfdx(diffRng,algRng)'*ldt);
lambda(end,diffRng)=ldt;
lambda(end,algRng)=lat;

%Alpha == 1 ==> Fully Implicit
alpha=1;
dfdx_tpdt=dfdx;
for ts=length(t)-1:-1:1
    dt=t(ts+1)-t(ts);
    %First solve the differential part
    dfdx_t=compute_dFdx(y(ts,:),p,rxParams);
    rhs=ldt-(1-alpha)*dt*...
        (dfdx_tpdt(diffRng,diffRng)'*ldt+dfdx_tpdt(algRng,diffRng)'*lat);
    A=eye(nDiff)+alpha*dt*...
        (dfdx_t(diffRng,diffRng)'-dfdx_t(algRng,diffRng)'*...
            (dfdx_t(algRng,algRng)'\dfdx_t(diffRng,algRng)'));
    ldt=A\rhs;
    lat=(dfdx_t(algRng,algRng)')\(-dfdx_t(diffRng,algRng)'*ldt);
    lambda(ts,diffRng)=ldt;
    lambda(ts,algRng)=lat;
    dfdx_tpdt=dfdx_t;
end