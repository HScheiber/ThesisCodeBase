function [TableFile_Out,C6,Energygrptables] = MakeTables(Settings,varargin)
    
    if isfield(Settings,'JobName') && ~isempty(Settings.JobName)
        TabName = [Settings.JobName '_Table'];
    else
        TabName = [Settings.Salt '_' Settings.Theory '_Table'];
    end
    
    p = inputParser;
    p.FunctionName = 'MakeTables';
    addOptional(p,'MDP_Minimize',false)
    addOptional(p,'TableName',TabName)
    addOptional(p,'SaveTables',true)
    parse(p,varargin{:});
    MDP_Minimize = p.Results.MDP_Minimize;
    TableName = p.Results.TableName;
    SaveTables = p.Results.SaveTables;
    
    [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
    
    switch Settings.Theory
        case {'JC' 'JC3P' 'JC4P' 'JCSD'}
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
            
            [U,C6] = JC_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'TF'
            [U,C6] = TF_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'HS'
            C6.MX = 1;C6.MM = 1;C6.XX = 1;
            U = HS_Potential_Generator(Settings,...
                'MDP_Minimize',true);
        case 'BH'
            [U,C6] = BH_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BD'
            [U,C6] = BD_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BE'
            [U,C6] = BE_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BF'
            [U,C6] = BF_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'Mie'
            [U,C6] = Mie_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        otherwise
            error(['Warning: Unknown model type: "' Settings.Theory '.'])
    end
    
    U_zero = zeros(size(U.r));
    ints = {'MX' 'MM' 'XX'};
    if Settings.Polarization
        subints = {'cc' 'sc' 'cs' 'ss'};
    else
        subints = {'cc'};
    end
    Energygrptables = cell(length(ints),length(subints));
    
    if SaveTables
        % First make blank table
        Uo = [U.r ; U_zero ; U_zero ; U_zero ; U_zero ; U_zero ; U_zero];
        U_out = deblank( sprintf(['%16.12e  %16.12e  %16.12e  %16.12e  %16.12e  %16.12e  %16.12e' newline],Uo(:)) );
        
        % Save table into current directory
        TableFile_Out = fullfile(Settings.WorkDir,[TableName '.xvg']);
        fid = fopen(TableFile_Out,'wt');
        fwrite(fid,regexprep(U_out,'\r',''));
        fclose(fid);
        
        for idx = 1:length(ints)
            int = ints{idx};
            for jdx = 1:length(subints)
                subint = subints{jdx};

                if strcmp(subint(1),'c') && strcmp(subint(2),'c') % core-core interactions
                    Uo = [U.r ; U.(int).f0 ; U.(int).df0 ; U.(int).g ; U.(int).dg ; U.(int).h ; U.(int).dh];
                else % shell interactions
                    % Output into gromacs format
                    Uo = [U.r ; U.(int).f0 ; U.(int).df0 ; U_zero ; U_zero ; U_zero ; U_zero];
                end
                
                Energygrptables{idx,jdx} = replace([int(1) '_' subint(1) ' ' int(2) '_' subint(2)],{'M' 'X' '_c'},{Metal Halide ''});
                U_out = deblank( sprintf(['%16.12e  %16.12e  %16.12e  %16.12e  %16.12e  %16.12e  %16.12e' newline],Uo(:)) );
                TableFile = fullfile(Settings.WorkDir,[TableName '_' strrep(Energygrptables{idx,jdx},' ','_') '.xvg']);
                
                % Save tables into current directory
                fid = fopen(TableFile,'wt');
                fwrite(fid,regexprep(U_out,'\r',''));
                fclose(fid);
            end
        end
    end
    
    Energygrptables = reshape(Energygrptables,1,[]);
    Energygrptables = sort(Energygrptables(~cellfun('isempty',Energygrptables)));
end