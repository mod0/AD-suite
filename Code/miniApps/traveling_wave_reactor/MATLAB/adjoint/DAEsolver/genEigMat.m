function M=genEigMat(f)

nCell=(length(f)-1)/4;
R=150;
B2geo=(2.405/R)^2;
h=400/nCell;
getMicros;

M=zeros(nCell);

%Process cross sections
Sf=msf9*f(1:nCell);
Sa=msa9*f(1:nCell)+...
    msa8*f(nCell+1:2*nCell)+...
    msa0*f(2*nCell+1:3*nCell);
St=mst9*f(1:nCell)+...
    mst8*f(nCell+1:2*nCell)+...
    mst0*f(2*nCell+1:3*nCell);
D=500./(3*St);
for i=1:nCell-1
    Dhalf(i,1)=((h/D(i)+h/D(i+1))/(2*h))^(-1);
end

M(1,1)=2.2*Sf(1)-Sa(1)-D(1)*B2geo-Dhalf(1)/h^2-2*D(1)/(h^2+4*D(1)*h);
M(1,2)=Dhalf(1)/h^2;
M(nCell,nCell)=2.2*Sf(nCell)-Sa(nCell)-D(nCell)*B2geo...
    -Dhalf(nCell-1)/h^2-2*D(nCell)/(h^2+4*D(nCell)*h);
M(nCell,nCell-1)=Dhalf(nCell-1)/h^2;
for i=2:nCell-1
    M(i,i)=2.2*Sf(i)-Sa(i)-D(i)*B2geo-Dhalf(i-1)/h^2-Dhalf(i)/h^2;
    M(i,i-1)=Dhalf(i-1)/h^2;
    M(i,i+1)=Dhalf(i)/h^2;
end

kill=12;
