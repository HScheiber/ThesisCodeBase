function TableFile_MX = MakeTables(Settings,varargin)

    p = inputParser;
    p.FunctionName = 'MakeTables';
    addOptional(p,'MDP_Minimize',false)
    parse(p,varargin{:});
    MDP_Minimize = p.Results.MDP_Minimize; % Gathers components of energy at the end in addition to the total energy.

    [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
    
    switch Settings.Theory(1:2)
        case 'JC'
            switch Settings.Theory
                case 'JC'
                    Settings.WaterModel = 'SPC/E';
                case 'JC3P'
                    Settings.WaterModel = 'TIP3P';
                case 'JC4P'
                    Settings.WaterModel = 'TIP4PEW';
                case 'JCSD'
                    Settings.WaterModel = 'SD';
            end
            
            [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,...
                'ReturnAsStructure','true','MDP_Minimize',MDP_Minimize);
        case 'TF'
            [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,...
                'ReturnAsStructure','true','MDP_Minimize',MDP_Minimize);
        case 'BH'
            [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,...
                'ReturnAsStructure','true','MDP_Minimize',MDP_Minimize);
    end

    % Find function minimum
    %% Grab the valleys of the potential
    TableName = [Settings.JobName '_Table'];
    Us = [U_MX,U_MM,U_XX];
    for idx = 1:3
        U = Us(idx);
        peaks_idx = [false islocalmax(U.Total(2:end))];
        peak_r = U.r(peaks_idx);
        if numel(peak_r) > 1
            peak_r = peak_r(1);
        end
        
        inflex_idx = [false islocalmax(U.dTotal(2:end)) | islocalmin(U.dTotal(2:end))];
        inflex_r = U.r(inflex_idx);
        
        if isempty(peak_r) || isempty(inflex_r) % Ensure peak exists
            add_wall = false;
        else
            inflex_r(inflex_r < peak_r) = [];
            inflex_r = inflex_r(1);
            add_wall = true;
        end
        
        if add_wall
            % Generate a steep repulsion beyond the peak
            below_peak_idx = (U.r < inflex_r);
            r = U.r(below_peak_idx); % nm
            D = 1e-15; % prefactor
            
            fwall = D./(r.^12) - D./(inflex_r.^12);
            dfwall = 12*D./(r.^13); % Wall -derivative
            
            % Kill the attractive interaction beyond the peak
            g_at_valley = U.g(U.r == inflex_r);
            N_belowpeak = sum(below_peak_idx);
            U.g(below_peak_idx) = repmat(g_at_valley,1,N_belowpeak);
            U.dg(below_peak_idx) = zeros(1,N_belowpeak);
            
            % Remove infinity at 0
            fwall(1) = fwall(2);
            dfwall(1) = 0;
            
            % Add this repulsion to the repulsive part of the function
            U.h(U.r < inflex_r) = U.h(U.r < inflex_r) + fwall;
            U.dh(U.r < inflex_r) = U.dh(U.r < inflex_r) + dfwall;
        end
        
%         % Testing
%         nm_per_m = 1e+9; % nm per m
%         NA = 6.0221409e23; % Molecules per mole
%         e_c = 1.60217662e-19; % Elementary charge in Coulombs
%         epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
%         k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
%         
%         if idx == 1
%             U.Total = -k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         else
%             U.Total =  k_0*(e_c^2).*(Settings.S.Q^2).*U.f0 + U.h + U.g ;
%         end
%         hold on
%         plot(U.r.*10,U.Total)
%         scatter(inflex_r.*10,U.Total(U.r == inflex_r))
%         ylim([-1000 1000])
        
        % Output into gromacs format
        Uo = [U.r ; U.f0 ; U.df0 ; U.g ; U.dg ; U.h ; U.dh];
        U_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],Uo(:)) );
        if idx == 1
            TableFile = fullfile(Settings.WorkDir,[TableName '.xvg']);
            TableFile_MX = TableFile;
        elseif idx == 2
            TableFile = fullfile(Settings.WorkDir,[TableName '_' Metal '_' Metal '.xvg']);
        elseif idx == 3
            TableFile = fullfile(Settings.WorkDir,[TableName '_' Halide '_' Halide '.xvg']);
        end
        
        % Save tables into current directory
        fid = fopen(TableFile,'wt');
        fwrite(fid,regexprep(U_out,'\r',''));
        fclose(fid);
    end

end