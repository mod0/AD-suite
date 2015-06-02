function dfdx = compute_dFdx(y,p,rxParams)

nCell=rxParams.nCell;
h=rxParams.h;
dfdx=zeros(length(y));

%Top 3*nCell rows are nuclide evolution equations
for i=1:nCell
    %N_9 Equations
    dfdx(i,i)=y(3*nCell+i)*p(2);
    dfdx(i,i+nCell)=-y(3*nCell+i)*p(4);
    dfdx(i,i+3*nCell)=-(y(i+nCell)*p(4)-y(i)*p(2));
    %N_8 Equations
    dfdx(i+nCell,i)=-y(3*nCell+i)*p(8)*p(2);
    dfdx(i+nCell,i+nCell)=y(3*nCell+i)*p(4);
    dfdx(i+nCell,i+3*nCell)=-(p(8)*y(i)*p(2)-y(i+nCell)*p(4));
    %N_0 Equations
    dfdx(i+2*nCell,i)=-2*y(i+3*nCell)*p(1);
    dfdx(i+2*nCell,i+2*nCell)=y(i+3*nCell)*p(6);
    dfdx(i+2*nCell,i+3*nCell)=-(2*y(i)*p(1)-y(i+2*nCell)*p(6));
end

%Now we need derivatives of the eigenproblem
%First do the nuclide densities in order
R=rxParams.R;
B2geo=rxParams.B2geo;

%We will loop thru interior cells first
for cl=2:nCell-1
    dvec(1)=p(10)/(3*(y(cl-1)*p(3)+y(nCell+cl-1)*p(5)+y(2*nCell+cl-1)*p(7)));
    dvec(2)=p(10)/(3*(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7)));
    dvec(3)=p(10)/(3*(y(cl+1)*p(3)+y(nCell+cl+1)*p(5)+y(2*nCell+cl+1)*p(7)));
    %Interior Cell, derivative wrt N_9
    dDdN=-p(10)/3*p(3)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
    dfdx(3*nCell+cl,cl)=(p(9)*p(1)-p(2)-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    %Interior cell, derivative wrt N_8
    dDdN=-p(10)/3*p(5)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
    dfdx(3*nCell+cl,nCell+cl)=(-p(4)-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    %Interior cell, derivative wrt N_0
    dDdN=-p(10)/3*p(7)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
    dfdx(3*nCell+cl,2*nCell+cl)=(-p(6)-B2geo*dDdN)*y(3*nCell+cl)...
        -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
        -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
    dfdx(3*nCell+cl-1,2*nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
    dfdx(3*nCell+cl+1,2*nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
end
%Now do cell 1, first wrt N_9
clear dvec
cl=1;
dvec(2)=p(10)/(3*(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7)));
dvec(3)=p(10)/(3*(y(cl+1)*p(3)+y(nCell+cl+1)*p(5)+y(2*nCell+cl+1)*p(7)));
dDdN=-p(10)/3*p(3)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,cl)=(p(9)*p(1)-p(2)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%Now N_8
dDdN=-p(10)/3*p(5)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,nCell+cl)=(-p(4)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%Now N_0
dDdN=-p(10)/3*p(7)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,2*nCell+cl)=(-p(6)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl+1,2*nCell+cl)=2/h^2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl+1));
%And now cell nCell, first wrt N_9
clear dvec
cl=nCell;
dvec(1)=p(10)/(3*(y(cl-1)*p(3)+y(nCell+cl-1)*p(5)+y(2*nCell+cl-1)*p(7)));
dvec(2)=p(10)/(3*(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7)));
dDdN=-p(10)/3*p(3)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,cl)=(p(9)*p(1)-p(2)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
%Now N_8
dDdN=-p(10)/3*p(5)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,nCell+cl)=(-p(4)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));
%Now N_0
dDdN=-p(10)/3*p(7)/(y(cl)*p(3)+y(nCell+cl)*p(5)+y(2*nCell+cl)*p(7))^2;
dfdx(3*nCell+cl,2*nCell+cl)=(-p(6)-B2geo*dDdN)*y(3*nCell+cl)...
    -2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1))...
    -2*h^2/(h^2+4*dvec(2)*h)^2*dDdN*y(3*nCell+cl);
dfdx(3*nCell+cl-1,2*nCell+cl)=2/h^2*dvec(1)^2/(dvec(1)+dvec(2))^2*dDdN*(y(3*nCell+cl)-y(3*nCell+cl-1));

%Now derivative of the eigenvalue problem wrt phi and Sigma_a^{ext}
M=genEigMat(y,p,rxParams);
dfdx(3*nCell+1:4*nCell,3*nCell+1:4*nCell)=M-eye(nCell)*y(end)/10^10;
dfdx(3*nCell+1:4*nCell,end)=-y(3*nCell+1:4*nCell);

%And derivatives of the scaling equation wrt N_9 and phi
dfdx(end,1:nCell)=y(3*nCell+1:4*nCell)*h*pi*R^2*p(1);
dfdx(end,3*nCell+1:4*nCell)=y(1:nCell)*h*pi*R^2*p(1);
        
    