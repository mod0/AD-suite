function [dgdp_pred dgdp_act]=testAdjoint(metFun)
% p=[msf9 msa9 mst9 msa8 mst8 msa0 mst0 Gamma nu alphaD];
pref=[1000.4e-24 (1026-7.968)*1e-24 1026e-24 ...
            500.717e-24 600.09e-24 ...
            5e-24 10e-24 ...
            .1*(1-1000.4/(1026-7.968)) 2.2 500]';
ppert=ones(size(pref))*1e-3;
ppert(8)=1e-3;

nCell=40;
[t0 y0 p0 rxParams0]=oneD_Wave(nCell,60*60*24*365*5,1,pref);
[g0 g_p g_xd g_xa]=metFun(t0,y0,p0,rxParams0);
lambda0=doAdjoint(t0,y0,p0,rxParams0,metFun);
dgdp_pred=compute_dGdp(lambda0,t0,y0,p0,rxParams0,metFun);
for k=1:length(pref)
    p=pref;
    p(k)=p(k)*(1+ppert(k));
    [t1 y1 p1 rxParams1]=oneD_Wave(nCell,60*60*24*365*5,1,p);
    [g1 g_p g_xd g_xa]=metFun(t1,y1,p1,rxParams1);
    dgdp_act(k)=(g1-g0)/(p1(k)-p0(k));
end

max(abs(dgdp_pred-dgdp_act)./dgdp_act)

pstring={'\sigma_{f,9}','\sigma_{a,9}','\sigma_{t,9}',...
    '\sigma_{a,8}','\sigma_{t,8}','\sigma_{a,0}','\sigma_{t,0}',...
    '\Gamma','\nu','\alpha_D'};
endLine='\\\hline';
for i=1:length(pref)
    fprintf('$%s$ & %.3e & %.3e & %.3e & %.3e %s\n',...
        char(pstring(i)), pref(i),dgdp_act(i),dgdp_pred(i),...
        abs((dgdp_pred(i)-dgdp_act(i))/dgdp_act(i)),endLine);
end