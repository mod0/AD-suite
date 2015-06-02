function [val eVec]=getSaAndFlux(macro,bndMeans)

% %Verification Numbers
% macro(1,:)=2.2*.6;
% macro(2,:)=.8;
% macro(3,:)=1/3;
% B2g=(pi/(4+4/3))^2;
% SS_sigA=macro(1,1)-macro(2,1)-macro(3,1)*B2g;
% bndMeans(1,:)=macro(3,1);


nCell=size(macro,2);
A=zeros(nCell);
h=bndMeans(2,1);
B2geo=(2.405/150)^2;

A(1,1)=macro(1,1)*h-h*macro(2,1)-h*macro(3,1)*B2geo-...
    bndMeans(1,1)/h-2*macro(3,1)/(h+4*macro(3,1));
A(1,2)=bndMeans(1,1)/h;
A(nCell,nCell)=macro(1,end)*h-h*macro(2,end)-h*macro(3,end)*B2geo-...
    bndMeans(1,end)/h-2*macro(3,end)/(h+4*macro(3,end));
A(nCell,nCell-1)=bndMeans(1,end)/h;
for row=2:nCell-1
    A(row,row-1)=bndMeans(1,row-1)/h;
    A(row,row+1)=bndMeans(1,row)/h;
    A(row,row)=macro(1,row)*h-macro(2,row)*h-macro(3,row)*h*B2geo-...
        bndMeans(1,row-1)/h-bndMeans(1,row)/h;
end
A=sparse(A/h);

% [V D]=eig(A);
% [val Ind]=max(sum(D));
% eVec=V(:,Ind);
% kill=12;
options.tol=1e-20;
[eVec val]=eigs(A,1,'la',options);

% plot(eVec)
% 
% x=ones(nCell,1);
% err=10;
% kold=1.2;
% while abs(err)>1e-6
% %     xnew=P*x;
% %     xnew=A\xnew;
%     xnew=A*x;
%     knew=norm(xnew)/norm(x);
%     xnew=xnew/knew;
%     err=norm(x-xnew);
%     x=xnew;
%     kold=knew;
% end
% eVal=knew;
% eVec2=xnew;
kill=12;