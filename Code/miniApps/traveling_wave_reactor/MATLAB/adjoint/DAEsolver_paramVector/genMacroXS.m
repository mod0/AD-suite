function [Sf Sa D Dhalf]=genMacroXS(f,p,rxParams)

nCell=rxParams.nCell;
h=rxParams.h;
%Process cross sections
Sf=p(1)*f(1:nCell);
Sa=p(2)*f(1:nCell)+...
    p(4)*f(nCell+1:2*nCell)+...
    p(6)*f(2*nCell+1:3*nCell);
St=p(3)*f(1:nCell)+...
    p(5)*f(nCell+1:2*nCell)+...
    p(7)*f(2*nCell+1:3*nCell);
% D=p(10)./(3*St);
% for i=1:nCell-1
%     Dhalf(i,1)=((h/D(i)+h/D(i+1))/(2*h))^(-1);
% end

D=St;
for i=1:nCell
    D(i)=p(10)/(3*St(i));
end
Dhalf=D;
Dhalf(end)=[];
for i=1:nCell-1
    Dhalf(i)=((h/D(i)+h/D(i+1))/(2*h))^(-1);
end