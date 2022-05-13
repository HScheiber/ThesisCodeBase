% Restart_TD_Metal_Complex
% A script to submit metal-complex calculations

% Some calculation parameters
Memory = '60GB';
Account = 'def-orvig';
Hours = 48; % Max calculation time
TD_roots = 10; % Number of excited states to calculate
Theories = {'CAM-B3LYP' 'wB97XD'};

% 'Ga' 'Fe' 'In' 'Sc' 'Lu'
Metal_Ions = {'Sc'};

% Available:
%                        'L8_Sym'                           'L8_Asym'
%             'HL7a_Sym' 'L7a_Sym'              'HL7a_Asym' 'L7a_Asym'
%             'HL7b_Sym' 'L7b_Sym'              'HL7b_Asym' 'L7b_Asym' 
% 'H2L6a_Sym' 'HL6a_Sym' 'L6a_Sym' 'H2L6a_Asym' 'HL6a_Asym' 'L6a_Asym'
%             'HL6b_Sym' 'L6b_Sym'              'HL6b_Asym' 'L6b_Asym'
% 'H2L6c_Sym' 'HL6c_Sym' 'L6c_Sym' 'H2L6c_Asym' 'HL6c_Asym' 'L6c_Asym'
% 'H2L5a_Sym' 'HL5a_Sym' 'L5a_Sym' 'H2L5a_Asym' 'HL5a_Asym' 'L5a_Asym'
%
% Submitted:
% 'L8_Sym' 'L8_Asym' 'L7a_Sym' 'L7a_Asym' 'L7b_Sym' 'L7b_Asym' 'L6a_Sym' 'L6a_Asym' 'L6b_Sym' 'L6b_Asym' 'L6c_Sym' 'L6c_Asym'
Structures = {'L7a_Asym'};

Mol = struct;

if ispc
    Template_Dir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\submission\Gaussian_scripts';
    pdir = pwd;
    Cpus = 32;
else
    [~,hostname] = system('hostname');
    switch lower(hostname(1:3))
        case 'gra'
            pdir = '/home/scheiber/project/HBEDpa_study';
            Cpus = 32;
        case 'ced'
            pdir = '/home/scheiber/orvig_project/HBEDpa_study';
            Cpus = 48;
    end
    Template_Dir = '/home/scheiber/ThesisCodeBase/submission/Gaussian_scripts';
end

% Load the template files
Input_Template = fileread( fullfile(Template_Dir,'Gaussian_TD_Restart.template') );
Subm_Template = fileread( fullfile(Template_Dir,'Gaussian_subm.template') );

% Add the global parameters
Input_Template = strrep(Input_Template,'##MEMORY##',Memory);
Input_Template = strrep(Input_Template,'##CPUS##',num2str(Cpus));
Input_Template = strrep(Input_Template,'##STATES##',num2str(TD_roots));

Subm_Template = strrep(Subm_Template,'##CPUS##',num2str(Cpus));
Subm_Template = strrep(Subm_Template,'##ACCOUNT##',Account);
Subm_Template = strrep(Subm_Template,'##TIME##',num2str(Hours));

for M_idx = 1:length(Metal_Ions)
    Mol.Metal_Ion = Metal_Ions{M_idx};
    [Mol,Multiplicities] = Get_Complex_Configs(Mol);
    
    % Add in basis set stuff
    Input_Template_M = strrep(Input_Template,'##BASIS_ROUTE_TD##',Mol.Basis_Route_TD);
    Input_Template_M = strrep(Input_Template_M,['##BASIS_SET_TD##' newline],Mol.Basis_Set_TD);
    
    if strcmp(Mol.Metal_Ion,'Lu')
        Input_Template_M = strrep(Input_Template_M,'##FULL##','full,');
    else
        Input_Template_M = strrep(Input_Template_M,'##FULL##','');
    end
    
    for T_idx = 1:length(Theories)
        Theory = Theories{T_idx};
        Input_Template_T = strrep(Input_Template_M,'##THEORY##',Theory);
    
        for MP_idx = 1:length(Multiplicities)
            Multiplicity = Multiplicities(MP_idx);
            
            if Multiplicity == 2
                Spin_State = '_LowSpin';
                continue % Skip low spin calculations
            elseif Multiplicity == 6
                Spin_State = '_HighSpin';
            else
                Spin_State = '';
            end

            for S_idx = 1:length(Structures)
                Mol.Structure = Structures{S_idx};
                Job_Name = [Mol.Metal_Ion '_' Mol.Structure Spin_State '_' Theory];

                % Add in Charge
                Input_Template_S = strrep(Input_Template_T,'##JOBNAME##',Job_Name);
                Subm_Template_S = strrep(Subm_Template,'##JOBNAME##',Job_Name);

                % Save the files and submit
                cdir = [pdir filesep Mol.Metal_Ion '_complexes' filesep strrep(Job_Name,['_' Theory],'')];
                if ~exist(cdir, 'dir') && ispc
                    mkdir(cdir);
                elseif ~exist(cdir, 'dir')
                   error(['No previous calculation found for: ' Job_Name]);
                end
                
                % Copy the checkpoint file
                if ~ispc
                    old_checkpoint = fullfile(cdir,[strrep(Job_Name,['_' Theory],'') '.chk']);
                    new_checkpoint = fullfile(cdir,[Job_Name '.chk']);
                    copyfile(old_checkpoint,new_checkpoint);
                end

                inp_text_filename = fullfile(cdir,[Job_Name '.com']);
                fid = fopen(inp_text_filename,'w');
                fprintf(fid,'%s',Input_Template_S);
                fclose(fid);

                submit_text_filename = fullfile(cdir,[Job_Name '.subm']);
                fid = fopen(submit_text_filename,'w');
                fprintf(fid,'%s',Subm_Template_S);
                fclose(fid);

                cd(cdir)
                disp(['submitting ' Job_Name]);
                if ~ispc
                    system(['sbatch ' submit_text_filename]);
                end
            end
        end
    end
end
cd(pdir)
