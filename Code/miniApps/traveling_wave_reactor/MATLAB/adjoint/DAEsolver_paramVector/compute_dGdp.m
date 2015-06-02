function dGdp=compute_dGdp(lambda,t,y,p,rxParams,metFun)
nAlg=rxParams.nAlg;
nDiff=rxParams.nDiff;
algRng=nDiff+1:(nAlg+nDiff);

[g g_p g_xd g_xa]=metFun(t,y,p,rxParams);
dfdx=compute_dFdx(y(end,:),p,rxParams);
dfdp=compute_dFdp(y(end,:),p,rxParams);
dGdp=g_p-g_xa*(dfdx(algRng,algRng)\dfdp(algRng,:));
dfdp_t=compute_dFdp(y(1,:),p,rxParams);
for i=1:(length(t)-1)
    dfdp_tpdt=compute_dFdp(y(i+1,:),p,rxParams);
    dGdp=dGdp-(t(i+1)-t(i))*.5*(lambda(i+1,:)*dfdp_tpdt+lambda(i,:)*dfdp_t);
    dfdp_t=dfdp_tpdt;
end
    