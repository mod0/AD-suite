function [g g_p g_xd g_xa]=EOL_fertileInv(t,y,p,rxParams);

g=rxParams.h*rxParams.R^2*pi*sum(y(end,rxParams.nCell+1:2*rxParams.nCell));
g_p=zeros(1,length(p));
g_xd=zeros(1,rxParams.nDiff);
g_xd(rxParams.nCell+1:2*rxParams.nCell)=rxParams.h*pi*rxParams.R^2;
g_xa=zeros(1,rxParams.nAlg);