clear all;
pm=0.4;
options = optimset('MaxIter',10,'TolX',1e-3,'TolFun', 1e-3,'Display','iter','GradObj','on','MaxFunEvals',30,'LargeScale','on');
A=[1];
b=[1.1];
[xx,fval,exitflag] = fmincon(@sens_check,pm,A,b,[],[],[],[],[],options);
xx
