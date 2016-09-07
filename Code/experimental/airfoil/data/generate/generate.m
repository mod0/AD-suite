function generate(I, J, flow)
    if(nargin < 2)
        I = 400;
        J = 600;
    end
    
    if(nargin < 3)
        flow = 0;
    end

    % generate the grid.dat file
    naca0012('old', I, J);

    % write independent variable to file.
    n_dim = 1;
    x_ind = [];
    x_ind(1) = 3.0 * atan(1.0) / 45.0;
    
    fid = fopen('x.dat', 'wt');
    for i = 1:n_dim
        fprintf(fid, ' %d\n', x_ind(i));
    end
    fclose(fid);
    
    % write all the parameters to 
    % parameters file.
    p_dim = 11;
    param = zeros(p_dim, 1);
    param(1) = 1.4d0;                                               % gam
    param(2) = 0.4d0;                                               % gm1
    param(3) = 0.9d0;                                               % cfl
    param(4) = 0.05d0;                                              % eps
    param(5) = 0.4d0;                                               % mach
    param(6) = 1.0d0;                                               % p
    param(7) = 1.0d0;                                               % r
    param(8) = sqrt(param(1) * param(6) / param(7)) * param(5);     % u
    param(9) = param(6)/(param(7) * param(2)) + 0.5 * param(8)^2;   % e
    param(10) = 20000;                                              % iter
    param(11) = 100;                                                % lift iter

    fid = fopen('parameters1.dat', 'wt');
    for i = 1:p_dim
        fprintf(fid, ' %d\n', param(i));
    end
    fclose(fid);
    
    if(flow ~= 0)
        ncell = 3 * I * J + 1;
        
        % compute the flow
        q = zeros(4, ncell);
        
        for i = 1:ncell
            q(1,i) = param(7);
            q(2,i) = param(7)*param(8);
            q(3,i) = 0.d0;
            q(4,i) = param(7)*param(9);
        end    
               
        fid = fopen('flow.dat','wt');
        for i = 1:ncell
            fprintf(fid,' %d %d %d %d\n',q(1, i), q(2, i), q(3, i), q(4, i));   
        end
        fclose(fid);
    else
        % write an empty flow.dat file
        fid = fopen('flow.dat', 'wt');
        fclose(fid);
    end
    
    % write dependent variable to file.
    m_dim = int32(param(10)/param(11));
    y_dep = zeros(m_dim, 1); 
    
    fid = fopen('y.dat', 'wt');
    for i = 1:m_dim
        fprintf(fid, ' %d\n', y_dep(i));
    end
    fclose(fid);
end
