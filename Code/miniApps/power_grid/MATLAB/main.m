clear;
global phimax pmmax
pm=0.4;
phimax = 1.0
pmmax = pm
tf11=[];
tf22=[];
f11=[];
f22=[];
options = optimset('MaxIter',10,'TolX',1e-3,'TolFun', 1e-3,'Display','iter','GradObj','on','MaxFunEvals',30,'LargeScale','on');
A=[1];
b=[1.1];
[xx,fval,exitflag] = fmincon(@sens_check,pm,A,b,[],[],[],[],[],options);
xx
% during optimization
[~,~,tf11, f11, ~,~] = sens_check(pmmax);
% after optimization
[~,~,tf22,f22,~,~] = sens_check(xx);
pt = plot(tf11,f11(:,1),tf22, f22(:,1),tf11, ones(size(tf11,1),1))
pt(1).LineWidth = 2;
pt(2).LineWidth = 2; 
pt(3).LineWidth = 2;
%plot(ones(1001,1));