function Output = Create_Liquid_Solid_Interface(Settings)

    Output.StructureChange = false;
    Output.SolidMelted = false;
    Output.LiquidFroze = false;
    Output.LiquidAmorphous = false;
    Output.Aborted = false;

    % Number of moles of liquid
    Settings.nmol_liquid = round((Settings.N_total/2)*Settings.Liquid_Fraction); % number of molecules needed
    
    % Find the approximate experimental density based on the temperature and pressure
    if ~isfield(Settings,'Ref_Density')
        Settings.Ref_Density = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
    end
    
    % Run a short equilibrating NPT simulation of the solid to obtain the correct solid density
    % This will update the SuperCell file
    if Settings.Equilibrate_Solid > 0
        Output = Equilibrate_Solid(Settings);
        if Output.Aborted
            return
        end
    end
    
    if Settings.Equilibrate_Liquid && Settings.GenCluster
        % Update liquid density to match average equilibrated density
        Output = Equilibrate_Liquid(Settings);
        if Output.Aborted
            return
        end
        Settings.Ref_Density = Output.Ref_Density;
        Minimize_Liquid_Cluster(Settings)
    elseif Settings.GenCluster
        if Settings.Verbose
            disp('Warning: Skipping Liquid Equilibration.')
        end
        Minimize_Liquid_Cluster(Settings)
    elseif Settings.Equilibrate_Liquid % Interface with pre-equilibration
        Output = Minimize_Equilibrate_Liquid_Interface(Settings);
        if Output.Aborted
            return
        end
    else % Interface without pre-equilibration
        Minimize_Liquid_Interface(Settings)
    end
end