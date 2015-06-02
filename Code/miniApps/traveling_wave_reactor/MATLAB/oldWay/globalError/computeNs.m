function nVec=computeNs(nVec,phiBnds,simTime,theta)
% Inputs:
%     nVec: vector containing intial values of nuclides
%         [U238 Pu239 FP]
%     phiBnds: length two vector containing initial and final value of phi 
%         during the timescale simTime.  Assume linear between initial and
%         final
%     simTime: total time to advance nuclide concentration
%     theta: time-stepping parameter.
%         theta=1     =>  Fully Implicit
%         theta=.5    =>  Crank Nicholson
%         theta=0     =>  Fully Explicit

%Check theta
if ((theta>1) || (theta<0))
    disp('Parameter Theta not in allowed bound, 0 to 1');
end

%Get Microscopic Cross Section Data
getMicros;

%Simulate Linear behavior for phi
slope=(phiBnds(2)-phiBnds(1))/simTime;
intercept=phiBnds(1);

%For now hard-code time step
dt=min([60*60 simTime]);

stillGoing=1;
count=0;
A=zeros(length(nVec));
rhs=zeros(length(nVec),1);
while stillGoing
    phi_n=count*dt*slope+intercept;
    phi_np1=(count+1)*dt*slope+intercept;
    if (count+1)*dt>=simTime
        dt=simTime-count*dt;
        stillGoing=0;
        phi_np1=phiBnds(2);
    end
    count=count+1;
    A(1,1)=1+dt*theta*phi_np1*msa9;
    A(1,2)=-theta*dt*phi_np1*msa8;
    A(1,3)=0;
    A(2,1)=-theta*Gamma*dt*msa9*phi_np1;
    A(2,2)=1+theta*dt*phi_np1*msa8;
    A(2,3)=0;
    A(3,2)=0;
    A(3,1)=-2*theta*dt*phi_np1*msf9;
    A(3,3)=1+theta*dt*phi_np1*msa0;
    rhs(1)=(1-theta)*dt*phi_n*(msa8*nVec(2)-msa9*nVec(1))+nVec(1);
    rhs(2)=dt*(1-theta)*phi_n*(Gamma*msa9*nVec(1)-msa8*nVec(2))+nVec(2);
    rhs(3)=dt*(1-theta)*phi_n*(2*msf9*nVec(1)-msa0*nVec(3))+nVec(3);
    nVec=A\rhs;
end
if sum(sum(nVec<0))>0
    warning('Negative Densities!')
end
    