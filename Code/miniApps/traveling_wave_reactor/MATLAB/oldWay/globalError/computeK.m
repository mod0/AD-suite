function [K dKdy dKdphi]=computeK(nVec,phi,sax)
nCell=length(phi);
h=400/nCell;
getMicros;
getMacros;
K=zeros(3*nCell,1);
dKdy=zeros(3*nCell);
dKdphi=zeros(3*nCell,nCell);
for i=1:nCell
    K(i)=phi(i)*(msa8*nVec(2,i)-msa9*nVec(1,i));
    K(nCell+i)=phi(i)*(Gamma*msa9*nVec(1,i)-msa8*nVec(2,i));
    K(2*nCell+i)=phi(i)*(2*msf9*nVec(1,i)-msa0*nVec(3,i));
end
% for i=1:nCell
%     if i==1
%         dphidt=nSpeed*((macro(1,i)-macro(2,i)-macro(3,i)*B2geo-sax)*phi(i)-...
%             bndMeans(1,i)/h^2*(phi(i)-phi(i+1))-...
%             2*macro(3,i)/(h^2+4*macro(3,i)*h)*phi(i));
%     elseif i==nCell
%         dphidt=nSpeed*((macro(1,i)-macro(2,i)-macro(3,i)*B2geo-sax)*phi(i)-...
%             bndMeans(1,i-1)/h^2*(phi(i)-phi(i-1))-...
%             2*macro(3,i)/(h^2+4*macro(3,i)*h)*phi(i));
%     else
%         dphidt=nSpeed*((macro(1,i)-macro(2,i)-macro(3,i)*B2geo-sax)*phi(i)-...
%             bndMeans(1,i-1)/h^2*(phi(i)-phi(i-1))-...
%             bndMeans(1,i)/h^2*(phi(i)-phi(i+1)));
%     end
%     dphidt=0;
%     dK(i)=phi(i)*(msa8*K(nCell+i)-msa9*K(i))+...
%         dphidt*K(i)/phi(i);
%     dK(nCell+i)=phi(i)*(Gamma*msa9*K(i)-msa8*K(nCell+i))+...
%         dphidt*K(nCell+i)/phi(i);
%     dK(2*nCell+i)=phi(i)*(2*msf9*K(i)-msa0*K(2*nCell+i))+...
%         dphidt*K(2*nCell+i)/phi(i);
% end
for i=1:nCell
    dKdy(i,i)=-phi(i)*msa9;
    dKdy(i,nCell+i)=phi(i)*msa8;
    dKdy(nCell+i,i)=Gamma*phi(i)*msa9;
    dKdy(nCell+i,nCell+i)=-phi(i)*msa8;
    dKdy(2*nCell+i,i)=2*msf9*phi(i);
    dKdy(2*nCell+i,2*nCell+i)=-phi(i)*msa0;
    dKdphi(i,i)=msa8*nVec(2,i)-msa9*nVec(1,i);
    dKdphi(nCell+i,i)=Gamma*msa9*nVec(1,i)-msa8*nVec(2,i);
    dKdphi(2*nCell+i,i)=2*msf9*nVec(1,i)-msa0*nVec(3,i);
end
kill=12;