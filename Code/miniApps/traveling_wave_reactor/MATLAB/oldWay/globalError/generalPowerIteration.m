function [eVal eVec]=generalPowerIteration(Sf,Sa,D,Dhalf,h)

nCell=length(Sf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Verification Lines: set up critical system with extrapolated boundaries
% slabWidth=nCell*h;
% macro(2,:)=0.8;
% macro(3,:)=1/3;
% rqdNuSf=(pi/(slabWidth+4*macro(3,1)))^2*macro(3,1)+macro(2,1);
% rqdNuSf=2.2*.6;
% macro(1,:)=rqdNuSf;
% bndMeans(1,:)=macro(3,1);
% bndMeans(2,:)=h;
% Dhalf(:)=macro(3,1);
% D(:)=1/3;
% Sa(:)=.8;
% Sf(:)=rqdNuSf;
% Sax=Sf(1)-Sa(1)-D(1)*(pi/(4+4*D(1)))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P=diag(h*Sf);
L=zeros(nCell);
B2geo=(2.405/150)^2;

for i=2:nCell-1
    L(i,i)=Dhalf(i-1)/h+Dhalf(i)/h+h*(Sa(i)+D(i)*B2geo);
    L(i,i-1)=-Dhalf(i-1)/h;
    L(i,i+1)=-Dhalf(i)/h;
end
L(1,1)=Dhalf(1)/h+2*D(1)/(h+4*D(1))+h*(Sa(1)+D(1)*B2geo);
L(1,2)=-Dhalf(1)/h;
L(nCell,nCell)=Dhalf(end)/h+2*D(end)/(h+4*D(end))+h*(Sa(end)+D(end)*B2geo);
L(nCell,nCell-1)=-D(end)/h;

counter=0;
x=ones(nCell,1);
err=10;
kold=1.2;
normold=norm(x);
maxold=max(abs(x));
xnew=x;


while abs(err)>3e-6
    counter=counter+1;
    xnew=P*x;
    xnew=L\xnew;
    newnorm=norm(xnew);
    knew=norm(xnew)/norm(x);
    xnew=xnew/knew;
    err=max(abs(xnew-x)./xnew);
    [counter err]
%     err=abs(norm(xnew/max(xnew))-norm(x/max(x)));
    x=xnew;
    kold=knew;
end
eVal=knew;
eVec=xnew;