% Submit_ion
%% Calculation parameters
Cutoff = 800;
Eps_default = 1e-15;
Eps_pgf_orb = 1e-20;
Eps_SCF = 1e-6;
vdw_cutoff = 700; % Only applies when using vdW or rVV10 corrections
smoothing = true; % uses SPLINE3 for XC derivative method
finer_grid = false; % Use finer grid for the XC
Hours = 3; % number of hours to give calculation
Cores = 25; % number of CPU cores to use
Box_Size = 15;
Broden_SCF = false;
restart = true; % Checks if job complete, if not then finds restart file and restarts the job
Increase_Radial_Grid = 0; % Increase radial grid size by this amount

% HFX options (hybrids only)
hfx.eps_schwarz = 1e-9;
hfx.max_memory = 3600;
hfx.cutoff_radius = Box_Size/2;

% Relativistic options
Relativistics = [true false]; % Switch to turn on relativistic Hamltonian
Rel.Method = 'DKH';
Rel.DKH_Order = 3;
Rel.Potential = 'FULL';
Rel.Transformation = 'ATOM';
Rel.ZORA_Type = 'MP'; % only applies to ZORA calculations
Rel.Z_Cutoff = 1;

% Options that change the system itself
%Ions = {'Li' 'Br' 'I'}; % 'Li' 'F' 'Cl' 'Br' 'I' 'At'
XC_Funcs =  {'TMTPSS-rVV10L'};
%{'optB88-vdW' 'PBE' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol-rVV10' ...
%         'rev-vdW-DF2' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS-rVV10L'};
% {'PBE' 'optB88' 'optB88-vdW' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol' 'PBEsol-rVV10' 'optPBE' 'optPBE-vdW' ...
%          'B86r' 'rev-vdW-DF2' 'PW86r' 'vdW-DF2' 'BEEF-vdW' ...
%          'SG4' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS' 'TMTPSS-rVV10L'};

Basis_Experiment = 'Sapporo-QZP'; % 'Sapporo-TZP' 'Sapporo-QZP'

%% Grab input templates
if ispc
    template_loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\CP2K_scripts';
    ALL_Potentials_loc = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\CP2K_scripts\ALL_POTENTIALS';
else
    template_loc = '/home/scheiber/ThesisCodeBase/CP2K_scripts';
    ALL_Potentials_loc = '/home/scheiber/ThesisCodeBase/CP2K_scripts/ALL_POTENTIALS';
    cd(['/project/6001647/scheiber/CP2K_Calcs/' Basis_Experiment])
end
pdir = pwd;

inp_template_txt = fileread([template_loc filesep 'input_atom_template.inp']);
subtext = fileread([template_loc filesep 'submission_template.subm']);

% Add in some universal parameters
inp_template_txt = strrep(inp_template_txt,'##CUTOFF##',num2str(Cutoff));
inp_template_txt = strrep(inp_template_txt,'##EPSDEFAULT##',num2str(Eps_default));
inp_template_txt = strrep(inp_template_txt,'##EPSPGFORB##',num2str(Eps_pgf_orb));
inp_template_txt = strrep(inp_template_txt,'##EPSSCF##',num2str(Eps_SCF));
inp_template_txt = strrep(inp_template_txt,'##BOXL##',num2str(Box_Size));
inp_template_txt = strrep(inp_template_txt,'##ALL_POTENTIAL_FILE##',ALL_Potentials_loc);

if Broden_SCF
    Broden = ['		&MIXING' newline ...
            '			ALPHA 0.5' newline ...
            '			METHOD BROYDEN_MIXING' newline ...
            '			NBUFFER 4' newline ...
            '			BROY_W0 0.0001' newline ...
            '		&END MIXING'];
    inp_template_txt = strrep(inp_template_txt,'##BRODEN##',Broden);
else
    inp_template_txt = strrep(inp_template_txt,'##BRODEN##','');
end

for rdx = 1:length(Relativistics)
    Relativistic = Relativistics(rdx);

    % Relativistic parameters
    if Relativistic
        Reltxt = cp2k_relativistic(Rel); %#ok<*UNRCH>
        inp_template_txt_rel = strrep(inp_template_txt,'##RELATIVISTIC##',Reltxt);
        func_add = ['-' Rel.Method];
        Ions = {'Li' 'Br' 'I'}; % 'Li' 'F' 'Cl' 'Br' 'I'
    else
        inp_template_txt_rel = strrep(inp_template_txt,'##RELATIVISTIC##','');
        func_add = '';
        Ions = {'Na' 'Li' 'F' 'Cl' 'Br' 'I'}; % 'Li' 'F' 'Cl' 'Br' 'I'
    end

    for idx = 1:length(XC_Funcs)
        XC_Func = XC_Funcs{idx};
        [XC_Func_input,Hybrid] = cp2k_XC_Func(XC_Func,vdw_cutoff,hfx,smoothing,finer_grid);

        inp_text_func = strrep(inp_template_txt_rel,'##XCFUNC##',XC_Func_input);

        for jdx = 1:length(Ions)
            Ion = Ions{jdx};

            if strcmp(Ion,'Li') || strcmp(Ion,'Na')
                Charge = '1';
            else
                Charge = '-1';
            end
            
            switch Basis_Experiment
                case 'pob-TZVP'
                    if Relativistic
                        switch Ion
                            case 'Li'
                                Basis_Name = 'def2-TZVPD';
                            otherwise
                                Basis_Name = 'Sapporo-DKH3-TZP-diffuse';
                        end
                    else
                        switch Ion
                            case {'Li' 'Na' 'F' 'Cl' 'Br'}
                                Basis_Name = 'def2-TZVPD';
                            case 'I'
                                Basis_Name = 'Sapporo-TZP-diffuse';
                        end
                    end
                case 'Sapporo-TZP'
                    if Relativistic
                        switch Ion
                            case {'Li' 'Na'}
                                Basis_Name = 'def2-TZVPD';
                            otherwise
                                Basis_Name = 'Sapporo-DKH3-TZP-diffuse';
                        end
                    else
                        switch Ion
                            case {'Li' 'Na'}
                                Basis_Name = 'def2-TZVPD';
                            otherwise
                                Basis_Name = 'Sapporo-TZP-diffuse';
                        end
                    end
                case 'Sapporo-QZP'
                    if Relativistic
                        switch Ion
                            case {'Li' 'Na'}
                                Basis_Name = 'def2-TZVPD';
                            otherwise
                                Basis_Name = 'Sapporo-DKH3-QZP-diffuse';
                        end
                    else
                        switch Ion
                            case {'Li' 'Na'}
                                Basis_Name = 'def2-TZVPD';
                            otherwise
                                Basis_Name = 'Sapporo-QZP-diffuse';
                        end
                    end
                otherwise
                    error(['Unknown Experiment: ' Basis_Experiment])
            end
            
            Basis_Set = cp2k_basis_set(Basis_Name,Ion,Increase_Radial_Grid);

            inp_text_ion = strrep(inp_text_func,'##BASISSET##',Basis_Set);
            inp_text_ion = strrep(inp_text_ion,'##CHARGE##',Charge);
            inp_text_ion = strrep(inp_text_ion,'##ATOM##',Ion);

            cdir = [pdir filesep XC_Func func_add filesep Ion];
            if ~exist(cdir, 'dir')
               mkdir(cdir)
            end

            Job_name = [XC_Func func_add '_' Ion];

            if restart
                is_complete = cp2k_check_complete(cdir,false);
                if is_complete
                    disp([Job_name ': Calculation completed, skipping.'])
                    continue
                else
                    d = dir([cdir filesep '*-1.restart']);

                    if isempty(d)
                        disp([Job_name ': Unable to find restart file, starting new.'])

                        submit_text_complete = strrep(subtext,'##JOBNAME##',Job_name);
                        submit_text_complete = strrep(submit_text_complete,'##HOURS##',num2str(Hours));
                        submit_text_complete = strrep(submit_text_complete,'##CORES##',num2str(Cores));
                        
                        inp_text_ion = strrep(inp_text_ion,'##JOBNAME##',Job_name);

                        inp_text_filename = fullfile(cdir,[Job_name '.inp']);
                        fid = fopen(inp_text_filename,'w');
                        fprintf(fid,'%s',inp_text_ion);
                        fclose(fid);
                    else
                        disp([Job_name ': Restart file found. Restarting Job.'])
                        % If multiple restart files exist, choose most recent
                        if size(d,1) > 1
                            Times = zeros(1,size(d,1));
                            for y = 1:size(d,1)
                                Times(y) = datenum(d(y).date,...
                                    'dd-mmm-yyyy HH:MM:SS');
                            end
                            [~,Idx] = max(Times);
                            d = d(Idx);
                        end

                        submit_text_complete = strrep(subtext,'##JOBNAME##',d.name);
                        submit_text_complete = strrep(submit_text_complete,'.inp','');
                        submit_text_complete = strrep(submit_text_complete,'-1.restart.out','.out');
                        submit_text_complete = strrep(submit_text_complete,'##HOURS##',num2str(Hours));
                        submit_text_complete = strrep(submit_text_complete,'##CORES##',num2str(Cores));
                    end
                end
            else
                submit_text_complete = strrep(subtext,'##JOBNAME##',Job_name);
                submit_text_complete = strrep(submit_text_complete,'##HOURS##',num2str(Hours));
                submit_text_complete = strrep(submit_text_complete,'##CORES##',num2str(Cores));
                inp_text_ion = strrep(inp_text_ion,'##JOBNAME##',Job_name);

                inp_text_filename = fullfile(cdir,[Job_name '.inp']);
                fid = fopen(inp_text_filename,'w');
                fprintf(fid,'%s',inp_text_ion);
                fclose(fid);
            end
            
            submit_text_complete = strrep(submit_text_complete,'##MEMORY##','3800M');
            submit_text_complete = strrep(submit_text_complete,'##CORESTOT##',num2str(Cores));

            submit_text_filename = fullfile(cdir,[Job_name '.subm']);
            fid = fopen(submit_text_filename,'w');
            fprintf(fid,'%s',submit_text_complete);
            fclose(fid);

            cd(cdir)
            disp(['submitting ' Job_name]);
            if ~ispc
                system(['sbatch ' submit_text_filename]);
                pause(1)
            end
        end
    end
end

cd(pdir)