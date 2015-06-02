function [val eVec]=getSaAndFlux(y,p,rxParams)
%A function which computes the value of an external absorber cross-section
%which makes a 1D system critical.  The value of the cross-section and the
%flux shape is returned



% %Verification Numbers
% macro(1,:)=2.2*.6;
% macro(2,:)=.8;
% macro(3,:)=1/3;
% B2g=(pi/(4+4/3))^2;
% SS_sigA=macro(1,1)-macro(2,1)-macro(3,1)*B2g;
% bndMeans(1,:)=macro(3,1);


% nCell=rxParams.nCell;
% A=zeros(nCell);
% h=rxParams.h;
% B2geo=rxParams.B2geo;



% %Construct the diffusion matrix
% A(1,1)=macro(1,1)*h-h*macro(2,1)-h*macro(3,1)*B2geo-...
%     bndMeans(1,1)/h-2*macro(3,1)/(h+4*macro(3,1));
% A(1,2)=bndMeans(1,1)/h;
% A(nCell,nCell)=macro(1,end)*h-h*macro(2,end)-h*macro(3,end)*B2geo-...
%     bndMeans(1,end)/h-2*macro(3,end)/(h+4*macro(3,end));
% A(nCell,nCell-1)=bndMeans(1,end)/h;
% for row=2:nCell-1
%     A(row,row-1)=bndMeans(1,row-1)/h;
%     A(row,row+1)=bndMeans(1,row)/h;
%     A(row,row)=macro(1,row)*h-macro(2,row)*h-macro(3,row)*h*B2geo-...
%         bndMeans(1,row-1)/h-bndMeans(1,row)/h;
% end
% A=sparse(A/h);

A=sparse(genEigMat(y,p,rxParams));

%This call returns the largest algebraic eigenvalue of the matrix A and the
%corresponding eigenvector
[eVec val]=eigs(A,1,'la');