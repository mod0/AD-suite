function y=lambdaPrimeTest(t,f,sol,p)

thisSol=deval(sol,t);
thisTrunc=computeTruncErr(sol,t,2);
thatTrunc=computeTruncErr(sol,t,2);
%The two ways seem to agree
if max(abs((thisTrunc-thatTrunc)./thisTrunc)) > 1e-5
    if ((abs(thisTrunc(2)*thatTrunc(2)) > 0) && (thisTrunc(2) > 1e-16))
        fprintf('Way 1 not equaling way 2!\n')
    end
end

y=zeros(4,1);
y(1)=thisSol(2)*p*f(1);
y(2)=thisSol(1)*p*f(1)+f(2);
y(3)=f(1)*thisTrunc(1)+f(2)*thisTrunc(2);
y(4)=thisSol(1)*thisSol(2)*f(1);