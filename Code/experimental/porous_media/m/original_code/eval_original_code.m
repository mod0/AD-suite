function [  ] = eval_original_code(data_directory, results_directory)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   input variables

    [n_dim, x]     = allocate_independent_variables(  );
    [m_dim, y]     = allocate_dependent_variables(  );
    [p_dim, param] = allocate_parameter_variables(  );
    
    [p_dim, param] = initialize_parameter_variables(   p_dim, param, data_directory );
    [n_dim, x]     = initialize_independent_variables( n_dim, x,     data_directory );

    [m_dim, y]     = evaluate_original_code(n_dim, m_dim, p_dim, x, y, param);
    
    save_dependent_variables( m_dim, y );

    [ n_dim, x ]     = deallocate_independent_variables( n_dim, x );
    [ n_dim, y ]     = deallocate_dependent_variables(   m_dim, y );
    [ npdim, param ] = deallocate_parameter_variables(   p_dim, param );
    return
end

function [ n_dim, x ] = allocate_independent_variables( )
% User-Application specific
% ===========================
% Standard AD-Suite Interface
% ===========================
    n_dim = 6;
    x     = zeros(n_dim,1);
    return
end

function [ m_dim, y ] = allocate_dependent_variables( )
% User-Application specific
% ===========================
% Standard AD-Suite Interface
% ===========================
    m_dim = 1;
    y     = zeros(m_dim,1);
    return
end

function [ p_dim, param ] = allocate_parameter_variables()
% User-Application specific
% ===========================
% N/A
% Standard AD-Suite Interface
% ===========================
% N/A
    p_dim = 0;
    param = [];
    return
end

function [ n_dim, x ] = initialize_independent_variables(n_dim, x)
    global data_directory

    % User-Application specific
    % ===========================
    % N/A
    % Standard AD-Suite Interface
    % =========================== 
    file  = strcat(data_directory, '/x.nc');
    n_dim = ncread(file,'n_dim');
    x     = ncread(file,'x');
    assert(isequal(n_dim,length(x)));
    return
end

function [ p_dim, param ] = initialize_parameter_variables(p_dim, param)
    global data_directory

    % User-Application specific
    % ===========================
    
    file1 = strcat(data_directory, '/parameters1.nc');
    file2 = strcat(data_directory, '/parameters2.nc');
    file3 = strcat(data_directory, '/parameters3.nc');
    file4 = strcat(data_directory, '/parameters4.nc');
    Grid.Nx   = double(ncread(file1,'NX'));
    Grid.Ny   = double(ncread(file1,'NY'));
    Grid.Nz   = double(ncread(file1,'NZ'));
    Grid.N    = Grid.Nx*Grid.Ny*Grid.Nz;
    Grid.hx   = ncread(file1,'hX');
    Grid.hy   = ncread(file1,'hY');
    Grid.hz   = ncread(file1,'hZ');
    Grid.V    = ncread(file1,'V');
    Fluid.vw  = ncread(file1,'vw');
    Fluid.vo  = ncread(file1,'vo');
    Fluid.swc = ncread(file1,'swc');
    Fluid.sor = ncread(file1,'sor');
    Steps.St  = double(ncread(file1,'St'));
    Steps.Pt  = double(ncread(file1,'Pt'));
    Steps.ND  = double(ncread(file1,'ND'));
    Steps.IR  = ncread(file1,'ir_const');
    Grid.K    = ncread(file2,'permeability');
    Por       = ncread(file3,'porosity');
    Grid.S    = ncread(file4,'saturation');

    Grid.por  = max(Por(:),1e-3);
    Grid.S    = Grid.S(:);
    Steps.IR  = 2*Steps.IR*Grid.V;
    % Standard AD-Suite Interface
    % ===========================
    param.Grid  = Grid;
    param.Fluid = Fluid;
    param.Steps = Steps;
    
    p_dim = size(param);
    return
end

function [ m_dim, y ] = evaluate_original_code(n_dim, m_dim, p_dim, x, y, param)
    addpath('./src')

    %Read Parameters
    Grid  = param.Grid;
    Fluid = param.Fluid;
    Steps = param.Steps;

    %Read independents
    mu    = x(1 : n_dim/2);
    sigma = x(n_dim/2+1 : end);

    %Run simulation
    Q = InitInOut(Grid, Steps.IR, mu, sigma);                    % Initialize Normal flow depending on x

    figure()
    Tt  = [];                                                    % For production curves
    Pc1 = []; 
    Pc2 = []; 
    res = [];

    S = Grid.S;
    for tp = 1:Steps.ND/Steps.Pt
        [~,V] = Pres(Grid,S,Fluid,Q);                              % Pressure solver
        for ts = 1:Steps.Pt/Steps.St;
            S         = NewtRaph(Grid, S, Fluid, V, Q, Steps.St);    % Implicit saturation solver
            Tt        = [Tt,(tp-1)*Steps.Pt+ts*Steps.St];            % Compute simulation time
            [Mw1,Mo1] = RelPerm(S(Grid.N), Fluid);                   % Compute output at well 1 and 2
            Mt1       = Mw1 + Mo1; 
            Pc1       = [Pc1,[Mw1/Mt1; Mo1/Mt1]];                    
            [Mw2,Mo2] = RelPerm(S(Grid.Ny),Fluid);             
            Mt2       = Mw2+Mo2; 
            Pc2       = [Pc2,[Mw2/Mt2; Mo2/Mt2]];
            res       = [res, -Mw1/Mt1+Mo1/Mt1+Mw2/Mt2-Mo2/Mt2];
            sum(Tt.*res)
            subplot(1,3,1)
            pcolor(reshape(S(Grid.Nx*Grid.Ny+1:2*Grid.Nx*Grid.Ny),Grid.Nx,Grid.Ny,1)');   % Plot saturation
            colorbar
            title('Visualization');
            caxis([0 1]) 
            shading flat
            sum((Tt(1:end)-[0 Tt(1:end-1)]) .*res)
            subplot(1,3,2)
            plot(Tt,Pc1(1,:))
            hold on 
            hold off
            title('Production at output 1');
            axis([0,Steps.ND,-0.05,1.05]);                 % Set correct axis
            legend('Water cut','Oil cut');                 % Set legend
            
            subplot(1,3,3)
            plot(Tt,Pc2(1,:))
            hold on 
            plot(Tt,Pc2(2,:));                         % Plot production data
            hold off
            title('Production at output 2');
            axis([0,Steps.ND,-0.05,1.05]);                 % Set correct axis
            legend('Water cut','Oil cut');                 % Set legend

            drawnow;                                           % Force update of plot
        end
    end

    % Write dependent Variables
    y(1) =  sum((Tt(1:end)-[0 Tt(1:end-1)]) .*res);
    return
end

function [ ] = save_dependent_variables(m_dim, y)
    global data_directory
    % Standard AD-Suite Interface
    % ===========================  
    % N/A
    % Standard AD-Suite Interface
    % ===========================
    file = strcat(data_directory, '/y.nc');  
    ncwrite(file, 'm_dim', m_dim);
    ncwrite(file,'y', y);
    return
end

function [ n_dim, x ] = deallocate_independent_variables( n_dim, x )
% User-Application specific
% ===========================
% N/A
% Standard AD-Suite Interface
% ===========================  
    n_dim = [];
    x     = [];
    return
end

function [ m_dim, y ] = deallocate_dependent_variables( m_dim, y )
% User-Application specific
% ===========================
% N/A
% Standard AD-Suite Interface
% ===========================  
    m_dim = [];
    y     = [];
    return
end

function [ p_dim, param  ] = deallocate_parameter_variables( p_dim, param )
% User-Application specific
% ===========================
% N/A
% Standard AD-Suite Interface
% ===========================
    p_dim = [];
    param = [];
    return
end