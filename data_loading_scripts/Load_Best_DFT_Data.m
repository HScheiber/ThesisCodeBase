function DFT = Load_Best_DFT_Data(varargin)

    % Parse optional inputs
    p = inputParser;
    p.FunctionName = 'Load_Best_DFT_Data';
    addOptional(p,'Conventional',true,@(x)validateattributes(x,...
        {'logical'},{'nonempty'}))
    addOptional(p,'Theory','TMTPSS-rVV10L',@(x)validateattributes(x,...
        {'char'},{'nonempty'}))
    parse(p,varargin{:})
    
    % Conversion units
    cm3_per_Ang3 = 1e-24; % cubic cm/cubic angstrom
    N_A = 6.02214076e23; % formula units/mol

    Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
    Molar_masses = [25.939 42.394 86.845 133.85 58.44];  % g/mol
    Structures = {'Rocksalt' 'Wurtzite' 'FiveFive' 'Sphalerite' 'NiAs' 'AntiNiAs' 'CsCl' 'BetaBeO'}; % 'BetaBeO'
    % Note: Best theory for LiI/CsCl is not corrected for scalar relativistic effects
    home = find_home;
    Data_Directory = [home filesep 'data'];
    CP2K_Data_Obj = load(fullfile(Data_Directory,'CP2K_Data.mat'));
    
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        MM = Molar_masses(idx);
        for jdx = 1:length(Structures)
            Structure = Structures{jdx};
            
            switch Salt
                case {'LiF' 'LiCl' 'NaCl'}
                    Theory = strrep(p.Results.Theory,'-','_');
                case 'LiI'
                    if strcmpi(Structure,'cscl')
                        Theory = strrep(p.Results.Theory,'-','_');
                    else
                        Theory = [strrep(p.Results.Theory,'-','_') '_DKH'];
                    end
                otherwise
                    Theory = [strrep(p.Results.Theory,'-','_') '_DKH'];
            end
            
            
            DFT.(Salt).(Structure).Energy = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).LE;
            DFT.(Salt).(Structure).a = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).a;
            DFT.(Salt).(Structure).b = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).b;
            DFT.(Salt).(Structure).c = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).c;
            
            if isnan(DFT.(Salt).(Structure).Energy)
                DFT.(Salt).(Structure).XC_E = nan;
                DFT.(Salt).(Structure).Non_Interacting_E = nan;
                DFT.(Salt).(Structure).Dispersion_E = nan;
                DFT.(Salt).(Structure).V = nan; % Volume per formula unit: Angstrom^3/Formula Unit
            else
                DFT.(Salt).(Structure).XC_E = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).XC_E; 
                DFT.(Salt).(Structure).Non_Interacting_E = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).Non_Interacting_E; 
                DFT.(Salt).(Structure).Dispersion_E = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).Dispersion_E;
                DFT.(Salt).(Structure).V = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).(Structure).Volume; % Angstrom^3/Formula Unit
            end
            
            DFT.(Salt).(Structure).density = MM/(DFT.(Salt).(Structure).V*cm3_per_Ang3*N_A); % g/cm^3
        end
        
        % If user selects primitive cell
        if ~p.Results.Conventional
            DFT.(Salt).Rocksalt.a = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).Rocksalt.a/sqrt(2);
            DFT.(Salt).Rocksalt.b = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).Rocksalt.b/sqrt(2);
            DFT.(Salt).Rocksalt.c = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).Rocksalt.c/sqrt(2);
            
            DFT.(Salt).FiveFive.a = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).FiveFive.c;
            DFT.(Salt).FiveFive.b = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).FiveFive.c;
            DFT.(Salt).FiveFive.c = CP2K_Data_Obj.Data.Sapporo_QZP.(Theory).(Salt).FiveFive.a;
        end
    end
end