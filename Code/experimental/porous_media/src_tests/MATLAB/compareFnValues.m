function compareFnValues()
    % First run runspe10.m in 2phase and save the output
    % Next run the fortran program and save the output
    % Now read the output from each
    load '../MATLAB/2phase/Pc.mat'
    PcMatlab = Pc';
    PcFortran = zeros(size(PcMatlab));
    PcFortran(:, 1) = readFortranOutput('../Fortran/DiscreteAdj/tempTAP/Pc1',401, 1);
    PcFortran(:, 2) = readFortranOutput('../Fortran/DiscreteAdj/tempTAP/Pc2',401, 1);
%     plot(PcFortran(:,1));
%     hold on;
%     plot(PcFortran(:,2));
%     plot(PcMatlab(:,1));
%     plot(PcMatlab(:,2));
    plot(PcMatlab - PcFortran);
    norm(PcMatlab - PcFortran);
%     legend('Fortran oil', 'Fortran water', 'Matlab oil', 'Matlab water');
end