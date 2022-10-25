function [TableFile_MX,C6] = MakeTables(Settings,varargin)

    p = inputParser;
    p.FunctionName = 'MakeTables';
    addOptional(p,'MDP_Minimize',false)
    addOptional(p,'TableName',Settings.JobName)
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
            
            [U_MX, U_MM, U_XX,C6] = JC_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'TF'
            [U_MX, U_MM, U_XX,C6] = TF_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'HS'
            C6.MX = 1;C6.MM = 1;C6.XX = 1;
            [U_MX, U_MM, U_XX] = HS_Potential_Generator(Settings,...
                'MDP_Minimize',true);
        case 'BH'
            [U_MX, U_MM, U_XX,C6] = BH_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BD'
            [U_MX, U_MM, U_XX,C6] = BD_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BE'
            [U_MX, U_MM, U_XX,C6] = BE_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'BF'
            [U_MX, U_MM, U_XX,C6] = BF_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        case 'Mie'
            [U_MX, U_MM, U_XX,C6] = Mie_Potential_Generator(Settings,...
                'MDP_Minimize',MDP_Minimize,...
                'Include_Dispersion_Scale',false);
        otherwise
            error(['Warning: Unknown model type: "' Settings.Theory '.'])
    end
    
	TableFile_MX = fullfile(Settings.WorkDir,[TableName '.xvg']);
    TableFile_MM = fullfile(Settings.WorkDir,[TableName '_' Metal '_' Metal '.xvg']);
    TableFile_XX = fullfile(Settings.WorkDir,[TableName '_' Halide '_' Halide '.xvg']);
    
    % Save tables into current directory
    if SaveTables
        fidMX = fopen(TableFile_MX,'wt');
        fwrite(fidMX,regexprep(U_MX,{'\r', '\n\n+'}',{'', '\n'}));
        fclose(fidMX);
        
        fidMM = fopen(TableFile_MM,'wt');
        fwrite(fidMM,regexprep(U_MM,{'\r', '\n\n+'}',{'', '\n'}));
        fclose(fidMM);
        
        fidXX = fopen(TableFile_XX,'wt');
        fwrite(fidXX,regexprep(U_XX,{'\r', '\n\n+'}',{'', '\n'}));
        fclose(fidXX);
    end

end