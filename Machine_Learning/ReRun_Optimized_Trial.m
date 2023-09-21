function ReRun_Optimized_Trial(Input_Settings)

    %% Load the model parameters
    if isstruct(Input_Settings)
        Settings = Input_Settings;
        batch_subm = false;
    elseif isfile(Input_Settings)
        Settings_dat = load(Input_Settings,'-mat');
        Settings = Settings_dat.Settings;
        batch_subm = true;
    else
        error('Input model does not exist, or is not a compatible data structure.')
    end


%% Final test of parameters on all structures for output
    % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
    if Settings.Parallel_Bayesopt || Settings.Parallel_LiX_Minimizer
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_Struct_Min = false;
        Settings.Parallel_LiX_Minimizer = true;
        Settings.MinMDP.Parallel_Min = false;
    elseif Settings.Parallel_Struct_Min
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_LiX_Minimizer = false;
        Settings.MinMDP.Parallel_Min = true;
    end
    Settings.Delete_Equil = false; % save the final MP calculation directories
    Settings.Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
    Settings.Verbose = true;
    Settings.CheckAmorphousHalide = true;
    [loss,~,UserData] = LiX_Minimizer(Settings,full_opt_point,...
        'Extra_Properties',true,'Therm_Prop_Override',true);
    
    ParNames = {params.Name};
    Pars = table;
    for idx = 1:numel(ParNames)
        Pars.(ParNames{idx}) = full_opt_point(idx);
    end
    
    % Coupled or each salt independent?
    [Salts,MultiMetal,MultiHalide] = Select_MultiSalt(Settings);
    
    % Load targets
    DFT = Load_Best_DFT_Data;
    Structures = Settings.Structures;
    N = length(Structures);
    if Settings.Loss_Options.Experimental_LE || Settings.Loss_Options.Experimental_LP
        Exp = Load_Experimental_Data;
        
        if Settings.Loss_Options.Experimental_LE
            for salt_idx = 1:numel(Salts)
                E_Correction = Exp.(Salts{salt_idx}).Rocksalt.E - DFT.(Salts{salt_idx}).Rocksalt.Energy;        
                for idx = 1:N
                    try
                        DFT.(Salts{salt_idx}).(Structures{idx}).Energy = DFT.(Salts{salt_idx}).(Structures{idx}).Energy + E_Correction;
                    catch
                        DFT.(Salts{salt_idx}).(Structures{idx}).Energy = nan;
                        DFT.(Salts{salt_idx}).(Structures{idx}).a = nan;
                        DFT.(Salts{salt_idx}).(Structures{idx}).b = nan;
                        DFT.(Salts{salt_idx}).(Structures{idx}).c = nan;
                    end
                end
            end
        end
        if Settings.Loss_Options.Experimental_LP
            for salt_idx = 1:numel(Salts)
                DFT.(Salts{salt_idx}).Rocksalt.a = Exp.(Salts{salt_idx}).Rocksalt.a_zero;
                DFT.(Salts{salt_idx}).Rocksalt.b = Exp.(Salts{salt_idx}).Rocksalt.b_zero;
                DFT.(Salts{salt_idx}).Rocksalt.c = Exp.(Salts{salt_idx}).Rocksalt.c_zero;
                DFT.(Salts{salt_idx}).Rocksalt.V = Exp.(Salts{salt_idx}).Rocksalt.V_zero;
                if isfield(Exp.(Salts{salt_idx}),'Wurtzite')
                    DFT.(Salts{salt_idx}).Wurtzite.a = Exp.(Salts{salt_idx}).Wurtzite.a_zero;
                    DFT.(Salts{salt_idx}).Wurtzite.b = Exp.(Salts{salt_idx}).Wurtzite.b_zero;
                    DFT.(Salts{salt_idx}).Wurtzite.c = Exp.(Salts{salt_idx}).Wurtzite.c_zero;
                    DFT.(Salts{salt_idx}).Wurtzite.V = Exp.(Salts{salt_idx}).Wurtzite.V_zero;
                end
            end
        end
    end
    
    for salt_idx = 1:numel(Salts)
        Minimization_Data = UserData.Minimization_Data;
        Finite_T_Data = UserData.Finite_T_Data;
        
        if MultiMetal || MultiHalide
            MinD_show = Minimization_Data.(Salts{salt_idx});
            FT_show = Finite_T_Data.(Salts{salt_idx});
        else
            MinD_show = Minimization_Data;
            FT_show = Finite_T_Data;
        end

        disp(repmat('*',1,80))
        disp(['Final Results - [Salt: ' Salts{salt_idx} '] - [Potential Form: ' Settings.Theory '] - [Model name: ' Settings.Trial_ID ']'])
        disp(repmat('*',1,80))

        format long g
        for idx = 1:length(Settings.Structures)
            disp(repmat('*',1,60))
            disp(Settings.Structures{idx})
            disp('Target (Exp/DFT) Data:')
            disp(repmat('*',1,60))
            disp(DFT.(Salts{salt_idx}).(Settings.Structures{idx}))
            disp(repmat('*',1,60))
            disp(Settings.Structures{idx})
            disp('Model Result:')
            disp(repmat('*',1,60))
            disp(MinD_show{idx})
        end
        disp(repmat('*',1,120))
        disp('Finite Temperature Data (Enthalpy in kJ/mol, Entropy in J/(mol K), Volume in A^3/Forumla Unit, Temperature in K)')
        disp(repmat('*',1,120))
        disp(FT_show)
        disp(repmat('*',1,120))

        disp(['Final Optimized Loss: ' num2str(loss,'%.10e')])
        disp('Final Optimized Parameters:')
        disp(Pars)
    end
    
    % Save final results
    save(Full_opt_filename,'full_opt_results','loss','full_opt_point',...
        'Minimization_Data','Finite_T_Data','Pars','Calculation_properties');
    
    %% Export final results to data store
    destination_folder = '/home/scheiber/project/Model_Building/Completed';
    inpfile = fullfile(Settings.OuterDir,[FileBase '.inp']);
    fullopt = fullfile(Settings.OuterDir,Full_opt_filename);
    bayesopt = fullfile(Settings.OuterDir,Results_filename);
    destfile = fullfile(destination_folder,Settings.Salt,[FileBase '_data.mat']);
    
    if isfile(fullopt)
        full_data = load(fullopt);
    else
        full_data = struct;
    end
    if isfile(bayesopt)
        bayesopt_dat = load(bayesopt);
    else
        bayesopt_dat = struct;
        bayesopt_dat.results = [];
    end
    
    if isfile(inpfile)
        input_model = load(inpfile,'-mat');
        full_data.Settings = input_model.Settings;
    end
    if isfile(Intermediate_Fullopt_file)
        fullopt_hist_dat = load(Intermediate_Fullopt_file);
        full_data.secondary_result = fullopt_hist_dat.intermediate_data;
    end
    full_data.bayesopt_results = bayesopt_dat.results;

    if ~isfolder(fullfile(destination_folder,Settings.Salt))
        mkdir(fullfile(destination_folder,Settings.Salt))
    end
    save(destfile,'full_data');
    diary off
    
end