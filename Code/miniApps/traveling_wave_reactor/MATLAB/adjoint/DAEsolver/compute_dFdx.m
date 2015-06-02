function dfdx = compute_dFdx(y)

getMicros;
nCell=(length(y)-1)/4;
h=400/nCell;
dfdx=zeros(length(y));

%Top 3*nCell rows are nuclide evolution equations
for i=1:nCell
    %N_9 Equations
    dfdx(i,i)=y(3*nCell+i)*msa9;
    dfdx(i,i+nCell)=-y(3*nCell+i)*msa8;
    dfdx(i,i+3*nCell)=-(y(i+nCell)*msa8-y(i)*msa9);
    %N_8 Equations
    dfdx(i+nCell,i)=-y(3*nCell+i)*Gamma*msa9;
    dfdx(i+nCell,i+nCell)=y(3*nCell+i)*msa8;
    dfdx(i+nCell,i+3*nCell)=-(Gamma*y(i)*msa9-y(i+nCell)*msa8);
    %N_0 Equations
    dfdx(i+2*nCell,i)=-2*y(i+3*nCell)*msf9;
    dfdx(i+2*nCell,i+2*nCell)=y(i+3*nCell)*msa0;
    dfdx(i+2*nCell,i+3*nCell)=-(2*y(i)*msf9-y(i+2*nCell)*msa0);
end

%Now we need derivatives of the eigenproblem
%First do the nuclide densities in order
nu=2.2;
R=150;
B2geo=(2.405/R)^2;
alphaD=500;
%We will loop thru interior cells first
for cl=2:nCell-1
    dvec(1)=alphaD/(3*(y(cl-1)*mst9+y(nCell+cl-1)*mst8+y(2*nCell+cl-1)*mst0));
    dvec(2)=alphaD/(3*(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0));
    dvec(3)=alphaD/(3*(y(cl+1)*mst9+y(nCell+cl+1)*mst8+y(2*nCell+cl+1)*mst0));
    %Interior Cell, derivative wrt N_9
    dDdN=-alphaD/3*mst9/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
    dfdx(3*nCell+cl,cl)=(nu*msf9-msa9-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    %Interior cell, derivative wrt N_8
    dDdN=-alphaD/3*mst8/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
    dfdx(3*nCell+cl,nCell+cl)=(-msa8-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    %Interior cell, derivative wrt N_0
    dDdN=-alphaD/3*mst0/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
    dfdx(3*nCell+cl,2*nCell+cl)=(-msa0-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,2*nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,2*nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
end
%Now do cell 1, first wrt N_9
clear dvec
cl=1;
dvec(2)=alphaD/(3*(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0));
dvec(3)=alphaD/(3*(y(cl+1)*mst9+y(nCell+cl+1)*mst8+y(2*nCell+cl+1)*mst0));
dDdN=-alphaD/3*mst9/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,cl)=(nu*msf9-msa9-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%Now N_8
dDdN=-alphaD/3*mst8/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,nCell+cl)=(-msa8-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%Now N_0
dDdN=-alphaD/3*mst0/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,2*nCell+cl)=(-msa0-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,2*nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%And now cell nCell, first wrt N_9
clear dvec
cl=nCell;
dvec(1)=alphaD/(3*(y(cl-1)*mst9+y(nCell+cl-1)*mst8+y(2*nCell+cl-1)*mst0));
dvec(2)=alphaD/(3*(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0));
dDdN=-alphaD/3*mst9/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,cl)=(nu*msf9-msa9-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
%Now N_8
dDdN=-alphaD/3*mst8/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,nCell+cl)=(-msa8-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
%Now N_0
dDdN=-alphaD/3*mst0/(y(cl)*mst9+y(nCell+cl)*mst8+y(2*nCell+cl)*mst0)^2;
dfdx(3*nCell+cl,2*nCell+cl)=(-msa0-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,2*nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));

%Now derivative of the eigenvalue problem wrt phi and Sigma_a^{ext}
M=genEigMat(y);
dfdx(3*nCell+1:4*nCell,3*nCell+1:4*nCell)=M-eye(nCell)*y(end)/10^10;
dfdx(3*nCell+1:4*nCell,end)=-y(3*nCell+1:4*nCell);

%And derivatives of the scaling equation wrt N_9 and phi
dfdx(end,1:nCell)=y(3*nCell+1:4*nCell)*h*pi*R^2*msf9;
dfdx(end,3*nCell+1:4*nCell)=y(1:nCell)*h*pi*R^2*msf9;
        
    