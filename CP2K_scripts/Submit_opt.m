% Submit_opt
%% Calculation parameters
Cutoff = 600;
kpoints = 6; % does not apply to hybrids
Eps_default = 1e-15;
Eps_pgf_orb = 1e-20;
Eps_SCF = 1e-5;
vdw_cutoff = 600; % Only applies when using vdW or rVV10 corrections
Ref_Theory = 'HSE06'; % for finding structure scale factor
smoothing = true; % uses SPLINE3 for XC derivative method
finer_grid = true; % Use finer grid for the XC
Hours = 48; % number of hours to give calculation
Cores = 25; % Number of CPU cores to use for calculation
Cores_total = 32; % Total number of CPU cores to request
Broden_SCF = false;
Fix_FC = true;
restart = true; % Checks if job complete, if not then finds restart file and restarts the job
skip_running = true; % checks if a job is already running and skips if it is.
Increase_Radial_Grid = 0; % Increase radial grid size by this amount
Fix_Ref_Cell = false;

% HFX options (hybrids only)
hfx.eps_schwarz = 1e-9;
hfx.max_memory = 3600;
hfx.cutoff_radius = 5;

% Relativistic options
Relativistics = [true false]; % Switch to turn on relativistic Hamltonian
Rel.Method = 'DKH';
Rel.DKH_Order = 3;
Rel.Potential = 'FULL';
Rel.Transformation = 'ATOM';
Rel.ZORA_Type = 'MP'; % only applies to ZORA calculations
Rel.Z_Cutoff = 1;

% Primary options that change the system and theory itself
XC_Funcs =  {'TMTPSS-rVV10L'};
% {'PBE' 'optB88' 'optB88-vdW' 'PBE-D3' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol' 'PBEsol-rVV10' 'optPBE' 'optPBE-vdW' ...
%         'B86r' 'rev-vdW-DF2' 'PW86r' 'vdW-DF2' 'BEEF-vdW' ...
%         'SG4' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS' 'TMTPSS-rVV10L'};
% {'optB88-vdW' 'PBE' 'PBE-D3' 'PBE-rVV10' 'PBE-rVV10L' 'PBEsol-rVV10' ...
%         'rev-vdW-DF2' 'SG4-rVV10m' 'SCAN' 'SCAN-rVV10' 'TMTPSS-rVV10L'};

Structures = {'BetaBeO'};
% {'Rocksalt' 'Wurtzite' 'FiveFive' 'NiAs' 'AntiNiAs' 'CsCl' 'Sphalerite'};
% 'BetaBeO'
Basis_Experiment = 'Sapporo-QZP'; % 'pob-TZVP' 'Sapporo-QZP' 'Sapporo-DZP'
cnt = 0;
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

inp_template_txt = fileread([template_loc filesep 'input_template.inp']);
subtext = fileread([template_loc filesep 'submission_template.subm']);

% Add in some universal parameters
inp_template_txt = strrep(inp_template_txt,'##CUTOFF##',num2str(Cutoff));
inp_template_txt = strrep(inp_template_txt,'##EPSDEFAULT##',num2str(Eps_default));
inp_template_txt = strrep(inp_template_txt,'##EPSPGFORB##',num2str(Eps_pgf_orb));
inp_template_txt = strrep(inp_template_txt,'##EPSSCF##',num2str(Eps_SCF));
inp_template_txt = strrep(inp_template_txt,'##ALL_POTENTIAL_FILE##',ALL_Potentials_loc);
inp_template_txt = strrep(inp_template_txt,'##RUNTYPE##','CELL_OPT');

if skip_running
    [~,jobstxt] = system('squeue -u $USER -o "%40j %20i %100Z" | grep -v "WORK_DIR" | grep "CP2K"');
    if ~isempty(jobstxt)
        current_job_info = textscan(strtrim(jobstxt),'%s %s %s','multipledelimsasone',true);
        current_job_names = current_job_info{1};
        current_job_ids = current_job_info{2};
        current_job_locs = current_job_info{3};
    else
        current_job_names = {};
        current_job_ids = {};
        current_job_locs = {};
    end
end

for rdx = 1:length(Relativistics)
    Relativistic = Relativistics(rdx);

    % Relativistic parameters
    if Relativistic
        Reltxt = cp2k_relativistic(Rel); %#ok<*UNRCH>
        inp_template_txt_rel = strrep(inp_template_txt,'##RELATIVISTIC##',Reltxt);
        func_add = ['-' Rel.Method];
        if strcmp(Basis_Experiment,'Sapporo-DZP')
            Salts = {'LiAt'};
        else
            Salts = {'LiBr' 'LiI'};
        end
    else
        inp_template_txt_rel = strrep(inp_template_txt,'##RELATIVISTIC##','');
        func_add = '';
        Salts = {'LiF' 'LiCl'}; % 'NaCl'
    end

    if Broden_SCF
        Broden = ['		&MIXING' newline ...
                '			ALPHA 0.5' newline ...
                '			METHOD BROYDEN_MIXING' newline ...
                '			NBUFFER 4' newline ...
                '			BROY_W0 0.0001' newline ...
                '		&END MIXING'];
        inp_template_txt_rel = strrep(inp_template_txt_rel,'##BRODEN##',Broden);
    else
        inp_template_txt_rel = strrep(inp_template_txt_rel,'##BRODEN##','');
    end

    for idx = 1:length(XC_Funcs)
        XC_Func = XC_Funcs{idx};
        [XC_Func_input,Hybrid,mGGA] = cp2k_XC_Func(XC_Func,vdw_cutoff,hfx,smoothing,finer_grid);

        inp_text_func = strrep(inp_template_txt_rel,'##XCFUNC##',XC_Func_input);
        if mGGA
            inp_text_func = strrep(inp_text_func,'##STRESSTENSOR##','DIAGONAL_NUMERICAL');
        else
            inp_text_func = strrep(inp_text_func,'##STRESSTENSOR##','ANALYTICAL');
        end

        for jdx = 1:length(Salts)
            Salt = Salts{jdx};
            [Metal,Halide] = Separate_Metal_Halide(Salt);
            
            switch Basis_Experiment
                case 'pob-TZVP'
                    Metal_Basis_Name = 'pob-TZVP'; % 'Sapporo-DKH3-TZP' 'Sapporo-TZP' 'pob-TZVP'
                    if Relativistic
                        Halide_Basis_Name = 'Sapporo-DKH3-TZP';
                    else
                        switch Halide
                            case {'F' 'Cl' 'Br'}
                                Halide_Basis_Name = 'pob-TZVP';
                            case 'I'
                                Halide_Basis_Name = 'Sapporo-TZP';
                        end
                    end
                case 'Sapporo-TZP'
                    Metal_Basis_Name = 'pob-TZVP'; % 'Sapporo-DKH3-TZP' 'Sapporo-TZP' 'pob-TZVP'
                    if Relativistic
                        Halide_Basis_Name = 'Sapporo-DKH3-TZP';
                    else
                        Halide_Basis_Name = 'Sapporo-TZP';
                    end
                case 'Sapporo-QZP'
                    Metal_Basis_Name = 'pob-TZVP';
                    if Relativistic
                        Halide_Basis_Name = 'Sapporo-DKH3-QZP';
                    else
                        Halide_Basis_Name = 'Sapporo-QZP';
                    end
                case 'Sapporo-DZP'
                    Metal_Basis_Name = 'pob-TZVP';
                    Halide_Basis_Name = 'Sapporo-DKH3-DZP-2012';
                otherwise
                    error(['Unknown Experiment: ' Basis_Experiment])
            end

            Basis_Set_M = cp2k_basis_set(Metal_Basis_Name,Metal,Increase_Radial_Grid);
            Basis_Set_X = cp2k_basis_set(Halide_Basis_Name,Halide,Increase_Radial_Grid);

            Basis_Set = [Basis_Set_M newline Basis_Set_X];

            inp_text_salt = strrep(inp_text_func,'##BASISSET##',Basis_Set);

            for kdx = 1:length(Structures)
                Structure = Structures{kdx};
                
                if strcmp(Salt,'LiAt') || strcmp(Structure,'BetaBeO')
                    Large_memory = true; % Use large memory nodes when true
                else
                    Large_memory = false; % Use large memory nodes when true
                end
                
                [Structure_Coord,ABC,Structure_Angle,Structure_Symmetry,Constraint] = ...
                    cp2k_Structures(Structure,Salt,Ref_Theory);

                [kpoints_txt,nrep] = cp2k_kpoints(kpoints,Hybrid,Structure,ABC,hfx.cutoff_radius);

                inp_text_structure = strrep(inp_text_salt,'##ANGLES##',Structure_Angle);
                inp_text_structure = strrep(inp_text_structure,'##COORDS##',Structure_Coord);
                inp_text_structure = strrep(inp_text_structure,'##SYMS##',Structure_Symmetry);
                inp_text_structure = strrep(inp_text_structure,'##ABCPAR##',ABC);
                inp_text_structure = strrep(inp_text_structure,'##KPOINTS##',kpoints_txt);
                inp_text_structure = strrep(inp_text_structure,'##REPA##',nrep.a);
                inp_text_structure = strrep(inp_text_structure,'##REPB##',nrep.b);
                inp_text_structure = strrep(inp_text_structure,'##REPC##',nrep.c);
                if Fix_FC
                    inp_text_structure = strrep(inp_text_structure,'##CONSTRAINT##',Constraint);
                else
                    inp_text_structure = strrep(inp_text_structure,['##CONSTRAINT##' newline],'');
                end

                cdir = [pdir filesep XC_Func func_add filesep Structure filesep Salt];
                if ~exist(cdir, 'dir')
                   mkdir(cdir)
                end

                Job_name = [XC_Func func_add '_' Structure '_' Salt];

                if skip_running
                    [job_running,jidx] = contained_in_cell(Job_name,current_job_names);
                    if job_running
                        job_running_loc = current_job_locs{jidx};
                        is_basis = ~isempty(regexp(job_running_loc,Basis_Experiment,'once'));
                    end
                    if job_running && is_basis
                        disp([Job_name ' (' Basis_Experiment '): Job already running with JOBID = ' current_job_ids{jidx} ', skipping job.'])
                        continue
                    end
                end

                if restart
                    is_complete = cp2k_check_complete(cdir,true);
                    if is_complete
                        disp([Job_name ': Calculation completed, skipping.'])
                        continue
                    else
                        d = dir([cdir filesep '*-1.restart']);
                        d_wf = dir([cdir filesep '*RESTART.kp']);

                        if isempty(d)
                            disp([Job_name ': Unable to find restart file, starting new.'])
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

                            restxt = fileread(fullfile(cdir,d.name));
                            abctxt = regexp(restxt,'&CELL\n(.+?)PERIODIC','tokens','once');
                            abcsc = textscan(abctxt{1},'%s %f %f %f');
                            ABC = [abcsc{1,2:end}];
                            a_vec = ABC(1,:);
                            b_vec = ABC(2,:);
                            c_vec = ABC(3,:);
                            a = num2str(norm(a_vec),'%10.15f');
                            b = num2str(norm(b_vec),'%10.15f');
                            c = num2str(norm(c_vec),'%10.15f');
                            FCtxt = regexp(restxt,'&COORD\n(.+?)SCALED','tokens','once');
                            FC = strtrim(FCtxt{1});

                            inp_text_structure = regexprep(inp_text_structure,'ABC +([0-9]|E|e|\.|-|+)+ +([0-9]|E|e|\.|-|+)+ +([0-9]|E|e|\.|-|+)+',...
                                ['ABC ' a ' ' b ' ' c],'once');
                            if ~Fix_FC
                                inp_text_structure = regexprep(inp_text_structure,'&COORD.+?&END COORD',...
                                    ['&COORD' newline '     SCALED' newline FC newline '    &END COORD'],'once');
                            end
                        end
                        
                        if isempty(d_wf)
                            disp([Job_name ': Unable to find WF restart file, starting new.'])
                        else
                            % If multiple restart files exist, choose most recent
                            if size(d_wf,1) > 1
                                Times = zeros(1,size(d_wf,1));
                                for y = 1:size(d_wf,1)
                                    Times(y) = datenum(d_wf(y).date,...
                                        'dd-mmm-yyyy HH:MM:SS');
                                end
                                [~,Idx] = max(Times);
                                d_wf = d_wf(Idx);
                            end

                            restxt = d_wf.name;

                            inp_text_structure = regexprep(inp_text_structure,'POTENTIAL_FILE_NAME (.+?)\n',...
                                ['POTENTIAL_FILE_NAME $1' newline '    WFN_RESTART_FILE_NAME ' restxt newline],'once');
                            inp_text_structure = regexprep(inp_text_structure,'SCF_GUESS +ATOMIC','SCF_GUESS RESTART','once');
                        end
                    end
                end

                
                if Fix_Ref_Cell
                    Ori_Cell = regexp(inp_text_structure,' *&CELL *\n.+?&END CELL','match','once');
                    
                    End_Cell = regexp(Ori_Cell,'\n *&END CELL','match','once');
                    End_Cell = End_Cell(2:end);
                    
                    Ref_Cell = strrep(Ori_Cell,'&CELL','&CELL_REF');
                    Ref_Cell = strrep(Ref_Cell,'&END CELL','&END CELL_REF');
                    Ref_Cell = regexprep(Ref_Cell,'(^|\n)    ','$1      ');
                    Ref_Cell = [Ref_Cell newline End_Cell];
                    New_Cell = strrep(Ori_Cell,End_Cell,Ref_Cell);
                    inp_text_structure = strrep(inp_text_structure,Ori_Cell,New_Cell);
                end
                
                submit_text_complete = strrep(subtext,'##JOBNAME##',Job_name);
                submit_text_complete = strrep(submit_text_complete,'##HOURS##',num2str(Hours));
                submit_text_complete = strrep(submit_text_complete,'##CORES##',num2str(Cores));
                submit_text_complete = strrep(submit_text_complete,'##CORESTOT##',num2str(Cores_total));
                if Large_memory
                    submit_text_complete = strrep(submit_text_complete,'##MEMORY##','16000M');
                else
                    submit_text_complete = strrep(submit_text_complete,'##MEMORY##','3800M');
                end
                
                inp_text_structure = strrep(inp_text_structure,'##JOBNAME##',Job_name);

                inp_text_filename = fullfile(cdir,[Job_name '.inp']);
                fid = fopen(inp_text_filename,'w');
                fprintf(fid,'%s',inp_text_structure);
                fclose(fid);

                submit_text_filename = fullfile(cdir,[Job_name '.subm']);
                fid = fopen(submit_text_filename,'w');
                fprintf(fid,'%s',submit_text_complete);
                fclose(fid);

                cd(cdir)
                disp(['submitting ' Job_name]);
                cnt = cnt+1;
                if ~ispc
                    system(['sbatch ' submit_text_filename]);
                    pause(1)
                end
            end
        end
    end
end

cd(pdir)