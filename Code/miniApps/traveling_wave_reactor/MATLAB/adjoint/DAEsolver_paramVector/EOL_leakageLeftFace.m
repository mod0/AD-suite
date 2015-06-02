function [gg g_p g_xd g_xa]=EOL_leakageLeftFace(t,y,p,rxParams);

nCell=rxParams.nCell;
params=gradientinit([y(end,:)';p]);
pRng=(4*nCell+1+1):length(params);

D1=params(pRng(10))/(3*(params(1)*params(pRng(3))+params(nCell+1)*params(pRng(5))+params(2*nCell+1)*params(pRng(7))));
g=2*D1*pi*rxParams.R^2*params(3*nCell+1)/(rxParams.h+4*D1);
gg=g.x;
g_p=full(g.dx(1,pRng));
g_xd=full(g.dx(1,1:rxParams.nDiff));
g_xa=full(g.dx(1,(rxParams.nDiff+1):(rxParams.nDiff+rxParams.nAlg)));
kill=12;