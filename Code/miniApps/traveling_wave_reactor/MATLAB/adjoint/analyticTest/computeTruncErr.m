function trunc=computeTruncErr(sol,x,way)

%for way==1, This code attempts to replicate exactly the way that the ode solver
%estimates local truncation error.  Here is a synopsis

    %If we are on the third or greater timestep, the code generates a cubic
    %interpolant between the most recent 4 computed solutions and
    %differentiates three times
     
    %If we are on the second time step, the code is using a BDF2 method, so
    %the local truncation error is proportional to the second derivative,
    %which is estimated by differencing the first derivative


if ((x<min(sol.x)) || (x>max(sol.x)))
    error('Cannot Extrapolate outside bounds of ODE Data');
end

if way==1
    topInd=find(x<=sol.x,1);
    if topInd==3
        t(1)=sol.x(1);
        y(:,1)=sol.y(:,1);
        t(3)=sol.x(2);
        y(:,3)=sol.y(:,2);
        t(2)=.5*(sol.x(1)+sol.x(2));
        y(:,2)=deval(sol,t(2));
        t(4)=x;
        y(:,4)=deval(sol,x);
        h=x-sol.x(2);
        trunc=doInterp(t,y,h);
        kill=12;
    elseif topInd==2  %The local truncation error is (1/2)y''*h^2
        %Compute this by differencing the scaled derivative
        [~, z]=deval(sol,sol.x(1));
        [~, znew]=deval(sol,x);
%         fprintf('Scaled Derivative?')
        h=x-sol.x(1);
        trunc=.5*(znew-z);
        kill=12;
    elseif topInd==1%at t==t_0 there is no truncation error -- set to 0 cubic
        trunc=zeros(size(sol.y,1),1);
    else
        t=sol.x([topInd-3:topInd]);
        t(end)=x;
        y=sol.y(:,[topInd-3:topInd]);
        thisy=deval(sol,x);
        y(:,end)=thisy;
        h=x-sol.x(topInd-1);
        trunc=doInterp(t,y,h);
        kill=12;
    end
        
%for way==2, I just use a cubic approximation as suggested in Mihai's notes
elseif way==2
    topInd=find(x<=sol.x,1);
    if topInd==1
        trunc=zeros(size(sol.y,1),1);
    else
        h=sol.x(topInd)-sol.x(topInd-1);
        myh=x-sol.x(topInd-1);
        if (myh>h)
            error('Someting Wong!');
        end
        err=sol.truncErr(topInd,:)';
        if topInd==2
            trunc=-1*err*(myh/h);
        else
            trunc=-1*err*(myh/h)^2;
        end
    end
end

end

function trunc=doInterp(t,y,h)
    c1=y(:,1)/((t(1)-t(2))*(t(1)-t(3))*(t(1)-t(4)));
    c2=y(:,2)/((t(2)-t(1))*(t(2)-t(3))*(t(2)-t(4)));
    c3=y(:,3)/((t(3)-t(1))*(t(3)-t(2))*(t(3)-t(4)));
    c4=y(:,4)/((t(4)-t(1))*(t(4)-t(2))*(t(4)-t(3)));
    est=(c1+c2+c3+c4)/2;
    trunc=est*abs(h)^2;
    p=pi;
    y0=1;
    tt=t(end);
    anld=(3*p^2*tt*y0)./(10000*exp((p*tt.^2)/200)) - (p^3*tt.^3*y0)./(1000000*exp((p*tt.^2)/200));
    anlEst=-1/12*anld*abs(h)^2;
    errEst=[trunc(1) anlEst];
    kill=12;
end