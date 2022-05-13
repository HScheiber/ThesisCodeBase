% Collect_Data

cd('/home/scheiber/project/CP2K_Calcs');
data_save_loc = '/home/scheiber/ThesisCodeBase/data/CP2K_Data.mat'; % Save the full dataset here
fileID = fopen('CP2K_Data_Summary.csv','wt');
Delta = 'Del ';
Epsilon = 'Er. ';
Ang = '(A)';
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
    'Theory',...
    'Salt',...
    [Epsilon 'LE[RS] Exp (%)'],...
    'LE[WZ]-LE[RS] (kJ/mol)',...
    'LE[5-5]-LE[RS] (kJ/mol)',...
    'LE[NiAs]-LE[RS] (kJ/mol)',...
    'LE[Spha]-LE[RS] (kJ/mol)',...
    [Epsilon 'a[RS] Exp (%)'],...
    [Epsilon 'a[WZ] Exp (%)'],...
    [Epsilon 'c[WZ] Exp (%)'],...
    [Epsilon 'c/a[WZ] Exp (%)'],...
    'MARE (%)'...
    );

Data = struct();
Warn_summary = '';
Missing_data = 0;

% List all basis experiment directories
disp('Beginning Search Through Raw Data...')
d_basis = dir('/home/scheiber/project/CP2K_Calcs');
dfolders_basis = d_basis([d_basis(:).isdir]);
Basis_dirs = dfolders_basis(~ismember({dfolders_basis(:).name},{'.','..'}));

for bdx = 1:length(Basis_dirs)
    Basis_name = Basis_dirs(bdx).name;
    Basis = strrep(Basis_name,'-','_');
    Basis_dir = [Basis_dirs(bdx).folder filesep Basis_name];
    disp('************************************')
    disp(['Current Basis Set Experiment: ' Basis_name])
    disp('************************************')
    Warn_summary = [Warn_summary newline '************************************' ...
        newline Basis_name newline '************************************']; %#ok<*AGROW>

    % List all theory directories
    d_theory = dir(Basis_dir);
    dfolders_theory = d_theory([d_theory(:).isdir]);
    XCFunc_dirs = dfolders_theory(~ismember({dfolders_theory(:).name},{'.','..'}));

    for idx = 1:length(XCFunc_dirs)
        Theory_name = XCFunc_dirs(idx).name;
        Theory = strrep(Theory_name,'-','_');
        XCFunc_dir = [XCFunc_dirs(idx).folder filesep Theory_name];
        disp(['Current Theory: ' Theory_name])

        % List all subdirectories
        d_struc = dir(XCFunc_dir);
        dfolders_struc = d_struc([d_struc(:).isdir]);
        Structure_Dirs = dfolders_struc(~ismember({dfolders_struc(:).name},{'.','..'}));

        for jdx = 1:length(Structure_Dirs)
            Structure = Structure_Dirs(jdx).name;
            Structure_dir = [XCFunc_dir filesep Structure];

            switch Structure_Dirs(jdx).name
                case {'Rocksalt' 'Wurtzite' 'FiveFive' 'NiAs' 'AntiNiAs' 'CsCl' 'BetaBeO' 'Sphalerite'} % Salts
                    % List all subdirectories
                    d_salt = dir(Structure_dir);
                    dfolders_salt = d_salt([d_salt(:).isdir]);
                    Salt_dirs = dfolders_salt(~ismember({dfolders_salt(:).name},{'.','..'}));

                    for kdx = 1:length(Salt_dirs)
                        Salt = Salt_dirs(kdx).name;
                        Salt_dir = [Structure_dir filesep Salt];

                        % Initialize data
                        Data.(Basis).(Theory).(Salt).(Structure).Energy = nan;
                        Data.(Basis).(Theory).(Salt).(Structure).a = nan;
                        Data.(Basis).(Theory).(Salt).(Structure).b = nan;
                        Data.(Basis).(Theory).(Salt).(Structure).c = nan;
                        Data.(Basis).(Theory).(Salt).(Structure).FC = '';

                        % Find all log files in current folder
                        d_out = dir([Salt_dir filesep '*.out']);

                        % Filter out slurm output files and grab the most
                        % recent log file only
                        if size(d_out,1) > 1
                            Times = zeros(1,size(d_out,1));
                            for y = 1:size(d_out,1)
                                if isempty(regexp(d_out(y).name,'slurm.*','ONCE'))

                                    Times(y) = datenum(d_out(y).date,'dd-mmm-yyyy HH:MM:SS');
                                end
                            end
                            [~,Idx] = max(Times);
                            d_out = d_out(Idx);
                        end
                        % Fix a naming issue
                        if ~isempty(regexp(d_out.name,'-1\.out','once'))
                            fixname = regexprep(d_out.name,'-1\.out','.out','once');
                            status = movefile(fullfile(d_out.folder,d_out.name),fullfile(d_out.folder,fixname),'f');
                            if status
                                d_out.name = fixname;
                            end
                        end

                        if isempty(d_out)
                            Warn_summary = [Warn_summary newline 'Missing log file for: ' Basis ' ' Theory ' ' Salt ' ' Structure];
                            Missing_data = Missing_data + 1;
                            continue
                        end

                        % Search log file for results of interest
                        jobtxt = fileread(fullfile(Salt_dir,d_out.name));
                        if contains(Basis_name,'fixed')
                            is_fixed = true;
                            comptxt = regexp(jobtxt,'SCF run converged in.+DBCSR STATISTICS','match','once');
                        else
                            is_fixed = false;
                            comptxt = regexp(jobtxt,'(GEOMETRY OPTIMIZATION COMPLETED|run CONVERGED!).+DBCSR STATISTICS','match','once');
                        end

                        if isempty(comptxt) && ~isfile(fullfile(Salt_dir,'CALCULATION_COMPLETE'))
                            
                            Warn_summary = [Warn_summary newline 'Calculation Not Complete For: ' Basis ' ' Theory ' ' Salt ' ' Structure];
                            Missing_data = Missing_data + 1;
                            
                            continue
%                             Entxt = regexp(jobtxt,'ENERGY\| Total FORCE_EVAL \( QS \) energy .a\.u\..: +(-|\.|[0-9]|E)+','tokens');
%                             Energies = cellfun(@str2double,[Entxt{:}]);
%                             TotE_cell = min(Energies);
                        else
                            if isfile(fullfile(Salt_dir,'CALCULATION_COMPLETE')) && isempty(comptxt)
                                comptxt = regexp(jobtxt,'SCF run converged in.+DBCSR STATISTICS','match','once');
                                is_fixed = true;
                            end
                            
                            Entxt = regexp(comptxt,'ENERGY\| Total FORCE_EVAL \( QS \) energy .a\.u\..: +(-|\.|[0-9]|E)+','tokens','once');
                            Overlap_Core_E = regexp(comptxt,'Overlap energy of the core charge distribution: +(-|\.|[0-9]|E)+','tokens','once');
                            Self_Core_E = regexp(comptxt,'Self energy of the core charge distribution: +(-|\.|[0-9]|E)+','tokens','once');
                            Core_Ham_E = regexp(comptxt,'Core Hamiltonian energy: +(-|\.|[0-9]|E)+','tokens','once');
                            Hartree_E = regexp(comptxt,'Hartree energy: +(-|\.|[0-9]|E)+','tokens','once');
                            XC_E = regexp(comptxt,'Exchange-correlation energy: +(-|\.|[0-9]|E)+','tokens','once');
                            Disp_E = regexp(comptxt,'Dispersion energy: +(-|\.|[0-9]|E)+','tokens','once');
                            GAPW_XC_E = regexp(comptxt,'GAPW\| Exc from hard and soft atomic rho1: +(-|\.|[0-9]|E)+','tokens','once');
                            GAPW_Core = regexp(comptxt,'GAPW\| local Eh = 1 center integrals: +(-|\.|[0-9]|E)+','tokens','once');
                            
                            TotE_cell = str2double(Entxt{1});
                            Overlap_Core_E = str2double(Overlap_Core_E{1});
                            Self_Core_E = str2double(Self_Core_E{1});
                            Core_Ham_E = str2double(Core_Ham_E{1});
                            Hartree_E = str2double(Hartree_E{1});
                            XC_E = str2double(XC_E{1});
                            if isempty(Disp_E)
                                Disp_E = 0;
                            else
                                Disp_E = str2double(Disp_E{1});
                            end
                            GAPW_XC_E = str2double(GAPW_XC_E{1});
                            GAPW_Core = str2double(GAPW_Core{1});
                        end
                        
                        if is_fixed
                            dres = dir([Salt_dir filesep '*.inp']);
                            if size(dres,1) > 1
                                Times = zeros(1,size(dres,1));
                                for y = 1:size(dres,1)
                                    Times(y) = datenum(dres(y).date,...
                                        'dd-mmm-yyyy HH:MM:SS');
                                end
                                [~,Idx] = max(Times);
                                dres = dres(Idx);
                            end
                            
                            if isempty(dres) % no restart file found
                                Warn_summary = [Warn_summary newline 'Missing restart file for: ' Basis ' ' Theory ' ' Salt ' ' Structure];
                                Missing_data = Missing_data + 1;
                                continue
                            end
                            
                            restxt = fileread(fullfile(Salt_dir,dres.name));
                            abctxt = regexp(restxt,'&CELL\n(.+?)PERIODIC','tokens','once');
                            abcsc = textscan(abctxt{1},'%s %f %f %f');
                            a_calc = abcsc{1,2}(1);
                            b_calc = abcsc{1,3}(1);
                            c_calc = abcsc{1,4}(1);
                            Input.alpha = abcsc{1,2}(2);
                            Input.beta = abcsc{1,3}(2);
                            Input.gamma = abcsc{1,4}(2);
                            alpha = Input.alpha;
                            beta = Input.beta;
                            gamma = Input.gamma;
                            TM = GenTransformMatrix(Input);
                            a_vec = TM(1,:).*a_calc;
                            b_vec = TM(2,:).*b_calc;
                            c_vec = TM(3,:).*c_calc;
                        else
                            dres = dir([Salt_dir filesep '*-1.restart']);
                            if size(dres,1) > 1
                                Times = zeros(1,size(dres,1));
                                for y = 1:size(dres,1)
                                    Times(y) = datenum(dres(y).date,...
                                        'dd-mmm-yyyy HH:MM:SS');
                                end
                                [~,Idx] = max(Times);
                                dres = dres(Idx);
                            end
                            if isempty(dres) % no restart file found
                                Warn_summary = [Warn_summary newline 'Missing restart file for: ' Basis ' ' Theory ' ' Salt ' ' Structure];
                                Missing_data = Missing_data + 1;
                                continue
                            end
                            
                            restxt = fileread(fullfile(Salt_dir,dres.name));
                            abctxt = regexp(restxt,'&CELL\n(.+?)PERIODIC','tokens','once');
                            abcsc = textscan(abctxt{1},'%s %f %f %f');
                            ABC = [abcsc{1,2:end}];
                            a_vec = ABC(1,:);
                            b_vec = ABC(2,:);
                            c_vec = ABC(3,:);
                            alpha = VecAngle(b_vec,c_vec);
                            beta = VecAngle(a_vec,c_vec);
                            gamma = VecAngle(a_vec,b_vec);
                            a_calc = norm(a_vec);
                            b_calc = norm(b_vec);
                            c_calc = norm(c_vec);
                        end

                        FCtxt = regexp(restxt,'&COORD\n(.+?)&END COORD','tokens','once');
                        FCtxt = strrep(FCtxt{1},'SCALED','');
                        FC = strrep(strtrim(FCtxt),newline,' ');

                        % Volume of unit cell in Angstrom^3
                        Volume_cell = (dot(cross(a_vec,b_vec),c_vec));

                        switch Structure
                            case 'Rocksalt'
                                if ismembertol(alpha,90,1e-3) % calculation performed in conventional cell
                                    a_conventional = a_calc;
                                    b_conventional = b_calc;
                                    c_conventional = c_calc;
                                    NF_per_Cell = 4;
                                else % In primitive cell
                                    a_conventional = a_calc*sqrt(2);
                                    b_conventional = b_calc*sqrt(2);
                                    c_conventional = c_calc*sqrt(2);
                                    NF_per_Cell = 1;
                                end
                            case 'Wurtzite'
                                a_conventional = a_calc;
                                b_conventional = b_calc;
                                c_conventional = c_calc;
                                NF_per_Cell = 2;
                            case 'CsCl'
                                a_conventional = a_calc;
                                b_conventional = b_calc;
                                c_conventional = c_calc;
                                NF_per_Cell = 1;
                            case 'FiveFive'
                                a_conventional = c_calc;
                                b_conventional = (3/sqrt(3))*a_calc;
                                c_conventional = a_calc;
                                a_primitive = a_calc;
                                b_primitive = b_calc;
                                c_primitive = c_calc;
                                NF_per_Cell = 2;
                            case 'BetaBeO'
                                a_conventional = a_calc;
                                b_conventional = b_calc;
                                c_conventional = c_calc;
                                NF_per_Cell = 4;
                            case 'Sphalerite'
                                if ismembertol(alpha,90,1e-3) % calculation performed in conventional cell
                                    a_conventional = a_calc;
                                    b_conventional = b_calc;
                                    c_conventional = c_calc;
                                    NF_per_Cell = 4;
                                else % In primitive cell
                                    a_conventional = a_calc*sqrt(2);
                                    b_conventional = b_calc*sqrt(2);
                                    c_conventional = c_calc*sqrt(2);
                                    NF_per_Cell = 1;
                                end
                            case {'NiAs' 'AntiNiAs'}
                                a_conventional = a_calc;
                                b_conventional = b_calc;
                                c_conventional = c_calc;
                                NF_per_Cell = 2;
                        end
                        Volume = Volume_cell/NF_per_Cell;               % Volume per ion pair
                        TotE = TotE_cell/NF_per_Cell;                   % Total Energy per ion pair
                                                                        % Energy Components per ion pair: %
                        Overlap_Core_E = Overlap_Core_E/NF_per_Cell;    % Overlap energy of the core charge distribution
                        Self_Core_E = Self_Core_E/NF_per_Cell;          % Self energy of the core charge distribution
                        Core_Ham_E = Core_Ham_E/NF_per_Cell;            % Core Hamiltonian energy
                        Hartree_E = Hartree_E/NF_per_Cell;              % Hartree energy
                        XC_E = XC_E/NF_per_Cell;                        % Exchange-correlation energy
                        Disp_E = Disp_E/NF_per_Cell;                    % Dispersion energy
                        GAPW_XC_E = GAPW_XC_E/NF_per_Cell;              % GAPW| Exc from hard and soft atomic rho1
                        GAPW_Core = GAPW_Core/NF_per_Cell;              % GAPW| local Eh = 1 center integrals
                        
                        Data.(Basis).(Theory).(Salt).(Structure).Energy = TotE;
                        Data.(Basis).(Theory).(Salt).(Structure).Volume = Volume;
                        Data.(Basis).(Theory).(Salt).(Structure).Overlap_Core_E = Overlap_Core_E;
                        Data.(Basis).(Theory).(Salt).(Structure).Self_Core_E = Self_Core_E;
                        Data.(Basis).(Theory).(Salt).(Structure).Core_Ham_E = Core_Ham_E;
                        Data.(Basis).(Theory).(Salt).(Structure).Hartree_E = Hartree_E;
                        Data.(Basis).(Theory).(Salt).(Structure).XC_E = XC_E + GAPW_XC_E;
                        Data.(Basis).(Theory).(Salt).(Structure).XC_Local_E = GAPW_XC_E;
                        Data.(Basis).(Theory).(Salt).(Structure).XC_Nonlocal_E = XC_E;
                        Data.(Basis).(Theory).(Salt).(Structure).Dispersion_E = Disp_E;
                        Data.(Basis).(Theory).(Salt).(Structure).GAPW_Core_E = GAPW_Core;
                        Data.(Basis).(Theory).(Salt).(Structure).Non_Interacting_E = Overlap_Core_E + ...
                            Self_Core_E + Core_Ham_E + Hartree_E + GAPW_Core;
                        Data.(Basis).(Theory).(Salt).(Structure).a = a_conventional;
                        Data.(Basis).(Theory).(Salt).(Structure).b = b_conventional;
                        Data.(Basis).(Theory).(Salt).(Structure).c = c_conventional;
                        Data.(Basis).(Theory).(Salt).(Structure).FC = FC;
                    end

                otherwise % Ions
                    Ion = Structure;

                    % Find all log files in current folder
                    d_out = dir([Structure_dir filesep '*.out']);

                    if isempty(d_out)
                        Warn_summary = [Warn_summary newline 'Unable to find log file for: ' Basis ' ' Theory ' ' Ion];
                        Missing_data = Missing_data + 1;
                        Data.(Basis).(Theory).(Ion).Energy = nan;
                        continue
                    end

                    TotE_cell = nan;
                    out_is_present = false;
                    % Search each log file for results of interest
                    for kdx = 1:length(d_out)

                        if ~isempty(regexp(d_out(kdx).name,'slurm.*','ONCE')) || ...
                                ~isempty(regexp(d_out(kdx).name,'8\.out','ONCE'))
                            continue
                        end
                        out_is_present = true;

                        jobtxt = fileread(fullfile(Structure_dir,d_out(kdx).name));
                        Entxt = regexp(jobtxt,'ENERGY\| Total FORCE_EVAL \( QS \) energy .a\.u\..: +(-|\.|[0-9]|E)+','tokens','once');
                        warntxt = regexp(jobtxt,'SCF run NOT converged','match','once');

                        if isempty(Entxt)
                            Warn_summary = [Warn_summary newline 'Unable to find ion energy for: ' Basis ' ' Theory ' ' Ion];
                            Missing_data = Missing_data + 1;
                            continue
                        elseif ~isempty(warntxt)
                            Warn_summary = [Warn_summary newline 'Ion energy did not satisfy SCF convergence for: ' Theory ' ' Ion];
                            Missing_data = Missing_data + 1;
                            continue
                        else
                            TotE_cell = str2double(Entxt{1});
                        end
                    end
                    if ~out_is_present
                        Warn_summary = [Warn_summary newline 'Unable to find log file for: ' Basis ' ' Theory ' ' Ion];
                        Missing_data = Missing_data + 1;
                        Data.(Basis).(Theory).(Ion).Energy = nan;
                        continue
                    end
                    Data.(Basis).(Theory).(Ion).Energy = TotE_cell; % Total energy of ion
            end
        end
    end
    Warn_summary = [Warn_summary newline '************************************'];
end
disp('Finished collecting Raw Data.')

% Experimental data
Experiment = Load_Experimental_Data;
E_Conv = 2625.4996394799; % Hartree -> kJ/mol conversion factor

%% Assigning lattice energy and Delta_E with experiment
disp('Calculating Lattice Energies and Experimental Differences.')
Basis_exps = fieldnames(Data);
for bdx = 1:length(Basis_exps)
    Basis = Basis_exps{bdx};
    disp(['Current Basis Experiment: ' strrep(Basis,'_','-')])

    % Calculate lattice energies, etc
    Theories = fieldnames(Data.(Basis));

    for idx = 1:length(Theories) % Loop through theories
        Theory = Theories{idx};
        disp(['Current Theory: ' strrep(Theory,'_','-')])
        Salts = fieldnames(Data.(Basis).(Theory));

        for jdx = 1:length(Salts) % Loop through theories
            Salt = Salts{jdx};
            [Metal,Halide] = Separate_Metal_Halide(Salt);

            if isempty(Metal) || isempty(Halide) % Skip over ion data
                continue
            end

            Structures = fieldnames(Data.(Basis).(Theory).(Salt));

            for kdx = 1:length(Structures)
                Structure = Structures{kdx};
                % Initialize data
                Data.(Basis).(Theory).(Salt).(Structure).LE = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_LE_Exp = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_LE_Exp_Percent = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_a_Exp = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_a_Exp_Percent = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_c_Exp = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_c_Exp_Percent = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_c_over_a_Exp = nan;
                Data.(Basis).(Theory).(Salt).(Structure).Delta_c_over_a_Exp_Percent = nan;

                % Lattice energy
                try
                    M_E = Data.(Basis).(strrep(Theory,'_D3','')).(Metal).Energy;
                    X_E = Data.(Basis).(strrep(Theory,'_D3','')).(Halide).Energy;
                catch
                    disp(['Metal or Halide energies unavailable for: ' Basis ' + ' strrep(Theory,'_D3','')]);
                    disp('Unable to calculate lattice energy')
                    M_E = nan;
                    X_E = nan;
                end
                Data.(Basis).(Theory).(Salt).(Structure).LE = ...
                    (Data.(Basis).(Theory).(Salt).(Structure).Energy - M_E - X_E)*E_Conv; % Lattice energy kj/mol

                % Comparison with experiment
                if strcmp(Structure,'Rocksalt') || strcmp(Structure,'Wurtzite')
                    % Error in LE compared with experiment (kJ/mol)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_LE_Exp = ...
                        (Data.(Basis).(Theory).(Salt).(Structure).LE - Experiment.(Salt).(Structure).E);

                    % Absolute error in LE compared with experiment (%)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_LE_Exp_Percent = ...
                        100*(Data.(Basis).(Theory).(Salt).(Structure).LE - ...
                        Experiment.(Salt).(Structure).E) / abs(Experiment.(Salt).(Structure).E);

                    % Error in lattice parameter [a] compared with experiment (Angstrom)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_a_Exp = ...
                        Data.(Basis).(Theory).(Salt).(Structure).a - ...
                        Experiment.(Salt).(Structure).a_zero;

                    % Absolute Error in lattice parameter [a] compared with experiment (%)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_a_Exp_Percent = ...
                        100*(Data.(Basis).(Theory).(Salt).(Structure).a - ...
                        Experiment.(Salt).(Structure).a_zero) / Experiment.(Salt).(Structure).a_zero;

                    % Error in lattice parameter [c] compared with experiment (Angstrom)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_c_Exp = ...
                        Data.(Basis).(Theory).(Salt).(Structure).c - ...
                        Experiment.(Salt).(Structure).c_zero;

                    % Absolute Error in lattice parameter [c] compared with experiment (%)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_c_Exp_Percent = ...
                        100*(Data.(Basis).(Theory).(Salt).(Structure).c - ...
                        Experiment.(Salt).(Structure).c_zero) / Experiment.(Salt).(Structure).c_zero;

                    % Error in c/a compared with experiment (unitless)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_c_over_a_Exp = ...
                        (Data.(Basis).(Theory).(Salt).(Structure).c/Data.(Basis).(Theory).(Salt).(Structure).a) - ...
                        Experiment.(Salt).(Structure).c_over_a;

                    % Absolute Error in c/a compared with experiment (%)
                    Data.(Basis).(Theory).(Salt).(Structure).Delta_c_over_a_Exp_Percent = ...
                        100*((Data.(Basis).(Theory).(Salt).(Structure).c/Data.(Basis).(Theory).(Salt).(Structure).a) - ...
                        Experiment.(Salt).(Structure).c_over_a) / Experiment.(Salt).(Structure).c_over_a;
                end

            end
            Data.(Basis).(Theory).(Salt).MARE = nanmean([ ...
                abs(Data.(Basis).(Theory).(Salt).Rocksalt.Delta_LE_Exp_Percent) ...
                abs(Data.(Basis).(Theory).(Salt).Rocksalt.Delta_a_Exp_Percent) ...
                %abs(Data.(Basis).(Theory).(Salt).Wurtzite.Delta_a_Exp_Percent) ...
                %abs(Data.(Basis).(Theory).(Salt).Wurtzite.Delta_c_Exp_Percent) ...
                ]);
                
        end
    end
end
disp('Finished Calculating Lattice Energies and Experimental Differences.')

% Assigning difference in lattice energy between structure and rocksalt at same level of theory
disp('Calculaing Differences in Lattice Energy Between Structures.')
for bdx = 1:length(Basis_exps)
    Basis = Basis_exps{bdx};
    disp(['Current Basis Experiment: ' strrep(Basis,'_','-')])
    Theories = fieldnames(Data.(Basis));

    for idx = 1:length(Theories) % Loop through theories
        Theory = Theories{idx};
        disp(['Current Theory: ' strrep(Theory,'_','-')])
        Salts = fieldnames(Data.(Basis).(Theory));

        for jdx = 1:length(Salts) % Loop through theories
            Salt = Salts{jdx};
            [Metal,Halide] = Separate_Metal_Halide(Salt);

            if isempty(Metal) || isempty(Halide) % Skip over ion data
                continue
            end

            Structures = fieldnames(Data.(Basis).(Theory).(Salt));

            for kdx = 1:length(Structures)
                Structure = Structures{kdx};
                
                if strcmp(Structure,'MARE') % Skip over MARE
                    continue
                end

                % Different in calculated LE and rocksalt calculated LE at same level of theory (kJ/mol)
                Data.(Basis).(Theory).(Salt).(Structure).Delta_LE_Rock = ...
                    Data.(Basis).(Theory).(Salt).(Structure).LE - Data.(Basis).(Theory).(Salt).Rocksalt.LE;

            end

            % Print out a data summary
            Out = cell(12,1);
            Out{1} = [strrep(Theory,'_','-') '/' strrep(Basis,'_','-')];
            Out{2} = Salt;
            try
                Out{3} = Data.(Basis).(Theory).(Salt).Rocksalt.Delta_LE_Exp_Percent;
            catch
                Out{3} = nan;
            end
            try
                Out{4} = Data.(Basis).(Theory).(Salt).Wurtzite.Delta_LE_Rock;
            catch
                Out{4} = nan;
            end
            try
                Out{5} = Data.(Basis).(Theory).(Salt).FiveFive.Delta_LE_Rock;
            catch
                Out{5} = nan;
            end
            try
                Out{6} = Data.(Basis).(Theory).(Salt).NiAs.Delta_LE_Rock;
            catch
                Out{6} = nan;
            end
            try
                Out{7} = Data.(Basis).(Theory).(Salt).Sphalerite.Delta_LE_Rock;
            catch
                Out{7} = nan;
            end
            try
                Out{8} = Data.(Basis).(Theory).(Salt).Rocksalt.Delta_a_Exp_Percent;
            catch
                Out{8} = nan;
            end
            try
                Out{9} = Data.(Basis).(Theory).(Salt).Wurtzite.Delta_a_Exp_Percent;
            catch
                Out{9} = nan;
            end
            try
                Out{10} = Data.(Basis).(Theory).(Salt).Wurtzite.Delta_c_Exp_Percent;
            catch
                Out{10} = nan;
            end
            try
                Out{11} = Data.(Basis).(Theory).(Salt).Wurtzite.Delta_c_over_a_Exp_Percent;
            catch
                Out{11} = nan;
            end
            try
                Out{12} = Data.(Basis).(Theory).(Salt).MARE;
            catch
                Out{12} = nan;
            end

            fprintf(fileID, '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',Out{:});

        end
    end
end
disp('Finished Calculaing Differences in Lattice Energy Between Structures.')
disp('Saving Data.')
Data.Experiment = Experiment;
fclose(fileID);

% Save the dataset
save(data_save_loc,'Data')
disp('Finished Saving Data! Job Complete.')
disp('************************************')
disp('Summary of missing data.')
disp('************************************')
disp(['Number of missing data points: ' num2str(Missing_data)])
disp(Warn_summary)
