function dfdp=compute_dFdp(y)

getMicros;
nCell=(length(y)-1)/4;
h=400/nCell;
dfdp=zeros(length(y),10);

% p=[msf9 msa9 mst9 msa8 mst8 msa0 mst0 Gamma nu alphaD]
%nuclide equations
for i=1:nCell
    %msf9
    dfdp(i,1)=0;
    dfdp(nCell+i,1)=0;
    dfdp(2*nCell+i,1)=-y(3*nCell+i)*2*y(i);
    %msa9
    dfdp(i,2)=y(3*nCell+i)*y(i);
    dfdp(nCell+i,2)=-y(3*nCell+i)*Gamma*y(i);
    dfdp(2*nCell+i,2)=0;
    %mst9
    dfdp(i,3)=0;
    dfdp(nCell+i,3)=0;
    dfdp(2*nCell+i,3)=0;
    %msa8
    dfdp(i,4)=-y(3*nCell+i)*y(nCell+i);
    dfdp(nCell+i,4)=y(3*nCell+i)*y(nCell+i);
    dfdp(2*nCell+i,4)=0;
    %mst8
    dfdp(i,5)=0;
    dfdp(nCell+i,5)=0;
    dfdp(2*nCell+i,5)=0;
    %msa0
    dfdp(i,6)=0;
    dfdp(nCell+i,6)=0;
    dfdp(2*nCell+i,6)=y(3*nCell+i)*y(2*nCell+i);
    %mst0
    dfdp(i,7)=0;
    dfdp(nCell+i,7)=0;
    dfdp(2*nCell+i,7)=0;
    %Gamma
    dfdp(i,8)=0;
    dfdp(nCell+i,8)=-y(3*nCell+i)*y(i)*msa9;
    dfdp(2*nCell+i,8)=0;
    %nu
    dfdp(i,9)=0;
    dfdp(nCell+i,9)=0;
    dfdp(2*nCell+i,9)=0;
    %alpha_D
    dfdp(i,10)=0;
    dfdp(nCell+i,10)=0;
    dfdp(2*nCell+i,10)=0;
end

%Now derivatives of the eigenvalue problem
nu=2.2;
R=150;
B2geo=(2.405/R)^2;
alphaD=500;
%msf9
dfdp(3*nCell+1:4*nCell,1)=y(3*nCell+1:4*nCell).*y(1:nCell)*nu;

%msa9
dfdp(3*nCell+1:4*nCell,2)=-y(3*nCell+1:4*nCell).*y(1:nCell);

%mst9
%Do the interior cells
for i=2:nCell-1
    dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
    dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
    dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
    dDids=-alphaD*y(i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
    dDimds=-alphaD*y(i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
    dDipds=-alphaD*y(i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
    dfdp(3*nCell+i,3)=-y(3*nCell+i)*B2geo*dDids...
        -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
        -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids);
end
%Now cell 1
i=1;
clear dvec;
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
dDids=-alphaD*y(i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDipds=-alphaD*y(i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
dfdp(3*nCell+1,3)=-y(3*nCell+1)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;
%Now cell nCell
i=nCell;
clear dvec;
dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dDids=-alphaD*y(i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDimds=-alphaD*y(i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
dfdp(3*nCell+i,3)=-y(3*nCell+i)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;

%msa8
dfdp(3*nCell+1:4*nCell,4)=-y(3*nCell+1:4*nCell).*y(nCell+1:2*nCell);

%mst8
%Do the interior cells
for i=2:nCell-1
    dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
    dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
    dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
    dDids=-alphaD*y(nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
    dDimds=-alphaD*y(nCell+i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
    dDipds=-alphaD*y(nCell+i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
    dfdp(3*nCell+i,5)=-y(3*nCell+i)*B2geo*dDids...
        -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
        -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids);
end
%Now cell 1
i=1;
clear dvec;
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
dDids=-alphaD*y(nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDipds=-alphaD*y(nCell+i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
dfdp(3*nCell+1,5)=-y(3*nCell+1)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;
%Now cell nCell
i=nCell;
clear dvec;
dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dDids=-alphaD*y(nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDimds=-alphaD*y(nCell+i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
dfdp(3*nCell+i,5)=-y(3*nCell+i)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;

%msa0
dfdp(3*nCell+1:4*nCell,6)=-y(3*nCell+1:4*nCell).*y(2*nCell+1:3*nCell);

%mst0
%Do the interior cells
for i=2:nCell-1
    dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
    dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
    dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
    dDids=-alphaD*y(2*nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
    dDimds=-alphaD*y(2*nCell+i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
    dDipds=-alphaD*y(2*nCell+i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
    dfdp(3*nCell+i,7)=-y(3*nCell+i)*B2geo*dDids...
        -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
        -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids);
end
%Now cell 1
i=1;
clear dvec;
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
dDids=-alphaD*y(2*nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDipds=-alphaD*y(2*nCell+i+1)/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0)^2);
dfdp(3*nCell+1,7)=-y(3*nCell+1)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipds+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;
%Now cell nCell
i=nCell;
clear dvec;
dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dDids=-alphaD*y(2*nCell+i)/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0)^2);
dDimds=-alphaD*y(2*nCell+i-1)/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0)^2);
dfdp(3*nCell+i,7)=-y(3*nCell+i)*B2geo*dDids...
    -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimds+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDids)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDids;

%Gamma
dfdp(3*nCell+1:4*nCell,8)=0;

%nu
dfdp(3*nCell+1:4*nCell,9)=y(3*nCell+1:4*nCell).*y(1:nCell)*msf9;

%alphaD
clear dDids dDimds dDdipds;
%Do the interior cells
for i=2:nCell-1
    dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
    dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
    dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
    dDidad=dvec(2)/alphaD;
    dDimdad=dvec(1)/alphaD;
    dDipdad=dvec(3)/alphaD;
    dfdp(3*nCell+i,10)=-y(3*nCell+i)*B2geo*dDidad...
        -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimdad+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDidad)...
        -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
            (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipdad+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDidad);
end
%Now cell 1
i=1;
clear dvec;
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dvec(3)=alphaD/(3*(y(i+1)*mst9+y(nCell+i+1)*mst8+y(2*nCell+i+1)*mst0));
dDidad=dvec(2)/alphaD;
dDipdad=dvec(3)/alphaD;
dfdp(3*nCell+1,10)=-y(3*nCell+1)*B2geo*dDidad...
    -(y(3*nCell+i)-y(3*nCell+i+1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(3))^2*dDipdad+2*dvec(3)^2/(dvec(2)+dvec(3))^2*dDidad)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDidad;
%Now cell nCell
i=nCell;
clear dvec;
dvec(1)=alphaD/(3*(y(i-1)*mst9+y(nCell+i-1)*mst8+y(2*nCell+i-1)*mst0));
dvec(2)=alphaD/(3*(y(i)*mst9+y(nCell+i)*mst8+y(2*nCell+i)*mst0));
dDidad=dvec(2)/alphaD;
dDimdad=dvec(1)/alphaD;
dfdp(3*nCell+i,10)=-y(3*nCell+i)*B2geo*dDidad...
    -(y(3*nCell+i)-y(3*nCell+i-1))/h^2*...
        (2*dvec(2)^2/(dvec(2)+dvec(1))^2*dDimdad+2*dvec(1)^2/(dvec(2)+dvec(1))^2*dDidad)...
    -2*y(3*nCell+i)*h^2/(h^2+4*dvec(2)*h)^2*dDidad;

%And finally derivatives of the power constraint -- only to msf9
dfdp(end,1)=sum(y(1:nCell).*y(3*nCell+1:4*nCell)*h*pi*R^2);
