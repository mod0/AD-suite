function B=compute_dMdN(cl,nuc,y)

nCell=(length(y)-1)/4;
nVec=zeros(3,nCell);
nVec(1,:)=y(1:nCell);
nVec(2,:)=y(nCell+1:2*nCell);
nVec(3,:)=y(2*nCell+1:3*nCell);
R=150;
B2geo=(2.405/R)^2;
h=400/nCell;
getMicros;
getMacros;

B=zeros(nCell);

switch nuc
    case 1
        sf=msf9;
        st=mst9;
        sa=msa9;
    case 2
        sf=0;
        st=mst8;
        sa=msa8;
    case 3
        sf=0;
        st=mst0;
        sa=msa0;
end
dDdN=-alphaD/3*st/(nVec(1,cl)*mst9+nVec(2,cl)*mst8+nVec(3,cl)*mst0)^2;
if ((cl>1) && (cl<nCell))
    B(cl,cl)=nu*sf-sa-B2geo*dDdN-...
        2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN-...
        2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
    B(cl,cl-1)=2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
    B(cl,cl+1)=2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
    B(cl-1,cl-1)=-2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
    B(cl-1,cl)=2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
    B(cl+1,cl+1)=-2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
    B(cl+1,cl)=2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
elseif cl==1
    B(cl,cl)=nu*sf-sa-B2geo*dDdN-...
        2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN-...
        2*h^2/(h^2+4*macro(3,cl)*h)^2*dDdN;
    B(cl,cl+1)=2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
    B(cl+1,cl)=2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
    B(cl+1,cl+1)=-2/h^2*macro(3,cl+1)^2/(sum(macro(3,cl:cl+1))^2)*dDdN;
elseif cl==nCell
    B(cl,cl)=nu*sf-sa-B2geo*dDdN-...
        2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN-...
        2*h^2/(h^2+4*macro(3,cl)*h)^2*dDdN;
    B(cl,cl-1)=2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
    B(cl-1,cl-1)=-2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
    B(cl-1,cl)=2/h^2*macro(3,cl-1)^2/(sum(macro(3,cl-1:cl))^2)*dDdN;
end