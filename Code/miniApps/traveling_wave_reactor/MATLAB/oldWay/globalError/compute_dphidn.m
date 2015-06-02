function [dphidn dMdN]=compute_dphidn(nVec)

nCell=size(nVec,2);
h=400/nCell;
getMicros;
getMacros;
dMdN=zeros(nCell,nCell,3*nCell);

B2geo=(2.405/150)^2;
dphidn=zeros(nCell,3*nCell);

A=genEigMat(nVec,B2geo);
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
% A=A/h;
[evecR evals]=eig(A);
[evecL evals]=eig(A.');evecL=conj(evecL);
for i=1:nCell
    evecR(:,i)=evecR(:,i)/norm(evecR(:,i));
    evecL(:,i)=evecL(:,i)/norm(evecL(:,i));
end
[val thisInd]=max(diag(evals));

nVec0=nVec;
macro0=macro;
sumVec=zeros(nCell,1);
eps=1e-8;
for nuc=1:3
    for cl=1:nCell
        nVec=nVec0;
        nVec(nuc,cl)=(1+eps)*nVec(nuc,cl);
        B=compute_dMdN(cl,nuc,nVec0,B2geo);
        dMdN(:,:,(nuc-1)*nCell+cl)=B;
        %Now compute expansion for u'_k (all except thisInd)
        %while at it, compute a_k
        ak=0;
        sumVec=zeros(nCell,1);
        for j=1:nCell
            if j==thisInd
                continue;
            end
            sj=evecL(:,j)'*evecR(:,j);
            aj=evecL(:,j)'*B*evecR(:,thisInd)/(sj*(evals(thisInd,thisInd)-evals(j,j)));
            ak=ak-aj*evecR(:,thisInd)'*evecR(:,j);
            sumVec=sumVec+aj*evecR(:,j);
        end
        dphidn(:,(nuc-1)*nCell+cl)=sumVec+ak*evecR(:,thisInd);
%         BB=(genEigMat(nVec,B2geo)-A)/(eps*nVec0(nuc,cl));
%         BBB=BB-B;
%         
%         if sum(abs(BBB(abs(BBB)>0)./B(abs(BBB)>0)))>.01
%             [B(abs(BBB)>0) BB(abs(BBB)>0) BBB(abs(BBB)>0) BBB(abs(BBB)>0)./B(abs(BBB)>0)]
%             [nuc cl]
%             kill=12;
%         end
%         kill=12;
    end
end
kill=12;
        