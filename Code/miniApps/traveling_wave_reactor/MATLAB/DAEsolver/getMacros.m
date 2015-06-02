%This script grabs the contents of nVec and makes cell-wise macroscopic
%cross sections.  The output is a matrix with 3 rows and nCell columns:
%     Row 1: nuSigmaF
%     Row 2: sigmaA
%     Row 3: D
%It also generates a matrix with 2 rows and (nCell-1) columns:
%     Row 1: Harmonic Mean of diffusion coefficient in cell i and i+1 (in column i)
%     Row 2: average spatial mesh between cell i and i+1 (in column i)
    
    
[nNuc nCell]=size(nVec);

nu=2.2;
macro(1,:)=nu*nVec(1,:)*msf9;
macro(2,:)=nVec(1,:)*msa9+nVec(2,:)*msa8+nVec(3,:)*msa0;
macro(3,:)=500./(3*(nVec(1,:)*mst9+nVec(2,:)*mst8+nVec(3,:)*mst0));


% %Verification tools
% macro(1,:)=.808;
% macro(2,:)=.8+sigmaX;
% macro(3,:)=1/3;

for i=1:nCell-1
    bndMeans(1,i)=((h/macro(3,i+1)+h/macro(3,i))/(2*h))^(-1);
    bndMeans(2,i)=h;
end

% [nNuc nCell]=size(nVec);
% 
% nu=2.2;
% macro(1,:)=nu*nVec(2,:)*msf9;
% macro(2,:)=nVec(1,:)*msa8+nVec(2,:)*msa9+nVec(3,:)*msa0+sigmaX;
% macro(3,:)=(.20)./(3*(nVec(1,:)*mst8+nVec(2,:)*mst9+nVec(3,:)*mst0));
% 
% 
% % %Verification tools
% % macro(1,:)=.808;
% % macro(2,:)=.8+sigmaX;
% % macro(3,:)=1/3;
% 
% for i=1:nCell-1
%     bndMeans(1,i)=((hVec(i+1)/macro(3,i+1)+hVec(i)/macro(3,i))/(hVec(i+1)+hVec(i)))^(-1);
%     bndMeans(2,i)=.5*(hVec(i+1)+hVec(i));
% end