function y=yprimeTest(t,f,p)
y=zeros(2,1);
y(1)=-p*f(2)*f(1);
y(2)=-t/200+f(2);