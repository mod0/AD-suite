function y=yprime(t,f,p,rxParams)
%This is the function which gives y'=f(y(t),t) for our DAE system
load tracker;
if t/(60*60*24*365)>tracker
    tracker=tracker+.5;
    fprintf('%.1f... ',t/(60*60*24*365));
    save tracker tracker;
end

nCell=rxParams.nCell;
y=zeros(length(f),1);

%Reactor parameters
L=rxParams.L;
R=rxParams.R;
B2geo=rxParams.B2geo;
h=rxParams.h;

%Target power is the targeted fission rate in the reactor or
%sum(phi*N_9*sigma_{f,9}*h)
targetPowerDensity=rxParams.targetPowerDensity;  %W/cm^3
joulesPerFission=rxParams.joulesPerFission;     %3.204e-11 = J/fission
%The target fission rate given the volume of the reactor
targetFissionRate=rxParams.targetFissionRate;


[Sf Sa D Dhalf]=genMacroXS(f,p,rxParams);
% %Process cross sections
% Sf=msf9*f(1:nCell);
% Sa=msa9*f(1:nCell)+...
%     msa8*f(nCell+1:2*nCell)+...
%     msa0*f(2*nCell+1:3*nCell);
% St=mst9*f(1:nCell)+...
%     mst8*f(nCell+1:2*nCell)+...
%     mst0*f(2*nCell+1:3*nCell);
% D=500./(3*St);
% for i=1:nCell-1
%     Dhalf(i,1)=((h/D(i)+h/D(i+1))/(2*h))^(-1);
% end

%Plutonium 239 densities
for i=1:nCell
    y(i)=f(i+3*nCell)*(p(4)*f(i+nCell)-p(2)*f(i));
end
%Uranium 238 densities
for i=nCell+1:2*nCell
    y(i)=f(i+2*nCell)*(p(8)*p(2)*f(i-nCell)-p(4)*f(i));
end
%Fission product densities
for i=2*nCell+1:3*nCell
    y(i)=f(i+nCell)*(2*p(1)*f(i-2*nCell)-p(6)*f(i));
end

%Algebraic Constraints
% %Flux vector
y(3*nCell+1)=(p(9)*Sf(1)*f(3*nCell+1)-...
        (f(4*nCell+1)/10^10+Sa(1)+D(1)*B2geo)*f(3*nCell+1)-...
        2*D(1)/(h^2+4*h*D(1))*f(3*nCell+1)-...
        Dhalf(1)/h^2*(f(3*nCell+1)-f(3*nCell+2)));
y(4*nCell)=(p(9)*Sf(nCell)*f(4*nCell)-...
        (f(4*nCell+1)/10^10+Sa(nCell)+D(nCell)*B2geo)*f(4*nCell)-...
        2*D(nCell)/(h^2+4*h*D(nCell))*f(4*nCell)-...
        Dhalf(nCell-1)/h^2*(f(4*nCell)-f(4*nCell-1)));
for i=2:nCell-1
    y(3*nCell+i)=(p(9)*Sf(i)*f(3*nCell+i)-...
        (f(4*nCell+1)/10^10+Sa(i)+D(i)*B2geo)*f(3*nCell+i)-...
        Dhalf(i-1)/h^2*(f(3*nCell+i)-f(3*nCell+i-1))-...
        Dhalf(i)/h^2*(f(3*nCell+i)-f(3*nCell+i+1)));
end

%Power Constraint
y(4*nCell+1)=sum(Sf.*f(3*nCell+1:4*nCell))*h*pi*R^2-targetFissionRate;
y(end);
kill=12;