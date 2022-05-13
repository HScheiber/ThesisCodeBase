function C6Damping = Init_C6Damping_Object

    % Dispersion damping
    C6Damping.MX = 0; % Set this to 0 for no dispersion damping!
    C6Damping.MM = 0;
    C6Damping.XX = 0;
    
    C6Damping.input_rvdw = false;
    C6Damping.rvdw.Li = 0.09; % nm
    C6Damping.rvdw.I = 0.20; % nm
    C6Damping.rvdw.Br = 0.18; % nm
    C6Damping.rvdw.Cl = 0.17; % nm
    C6Damping.rvdw.F = 0.12; % nm
    
    % Damping on additional function
    C6Damping.N.MX = 0; % Set this to 0 for no damping!
    C6Damping.N.MM = 0;
    C6Damping.N.XX = 0;
end