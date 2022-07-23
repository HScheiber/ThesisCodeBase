function Create_Liquid_Solid_Interface(Settings)

    % Number of moles of liquid
    Settings.nmol_liquid = round((Settings.N_total/2)*Settings.Liquid_Fraction); % number of molecules needed
    
    % Find the approximate experimental density based on the temperature and pressure
    if ~isfield(Settings,'Ref_Density')
        warning('off','MATLAB:UndefinedFunction')
        Settings.Ref_Density = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
        warning('on','MATLAB:UndefinedFunction')
    end
    
    % Run a short equilibrating NPT simulation of the solid to obtain the correct solid density
    % This will update the SuperCell file
    if Settings.Equilibrate_Solid > 0
        Equilibrate_Solid(Settings)
    end
    
    if Settings.Equilibrate_Liquid && Settings.GenCluster
        % Update liquid density to match average equilibrated density
        Settings.Ref_Density = Equilibrate_Liquid(Settings);
        Minimize_Liquid_Cluster(Settings)
    elseif Settings.GenCluster
        disp('Warning: Skipping Liquid Equilibration.')
        Minimize_Liquid_Cluster(Settings)
    elseif Settings.Equilibrate_Liquid % Interface with pre-equilibration
        Minimize_Equilibrate_Liquid_Interface(Settings);
    else % Interface without pre-equilibration
        Minimize_Liquid_Interface(Settings)
    end

end