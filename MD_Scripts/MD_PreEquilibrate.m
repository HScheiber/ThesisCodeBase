function Output = MD_PreEquilibrate(Directory)

Output.StructureChange = false;
Output.SolidMelted = false;
Output.LiquidFroze = false;
Output.LiquidAmorphous = false;
Output.Aborted = false;

Settings = load(fullfile(Directory,'TempJobInfo.mat'));

% Place "Settings" sub-structure into the main input structure
if isfield(Settings,'Settings')
    f = fieldnames(Settings.Settings);
    for i = 1:length(f)
        Settings.(f{i}) = Settings.Settings.(f{i});
    end
    Settings = rmfield(Settings,'Settings');
end
Inp_Settings = Settings;
Settings.WorkDir = Directory;
if ~isfield(Settings,'Verbose')
    Settings.Verbose = true;
end

if Settings.Verbose
    disp('*** Fast Pre-Equilibration Selected ***')
end
% Set the number of steps
timesteps = Settings.PreEquilibration/Settings.MDP.dt;
Compressibility = Settings.QECompressibility.*ones(1,length(Settings.Target_P)); % bar^(-1)
if length(Compressibility) > 2
    Compressibility(4:6) = 0;
end
Compressibility = regexprep(num2str(Compressibility),' +',' ');
tau_p = Settings.MDP.dt; % ps
tau_t = Settings.MDP.dt; % ps
    
nstpcouple = max(round(tau_p/(20*Settings.MDP.dt)),1);
nsttcouple = max(round(tau_t/(20*Settings.MDP.dt)),1);

Target_P = regexprep(num2str(Settings.Target_P),' +',' ');
Target_T = num2str(Settings.Target_T);

% Ensure fast equilibration with Berendsen barostat + small time constant
MDP_Template = regexprep(Settings.MDP_Template,'\r\n','\n');
MDP_Template = regexprep(MDP_Template,'(nsteps += *)(.+?)( *);',['$1' num2str(timesteps) '$3;']);
MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$1100$3;');
MDP_Template = regexprep(MDP_Template,'(annealing += *)(.+?)( +);','$1no$3;');

MDP_Template = regexprep(MDP_Template,'(\; Pressure Coupling\n).+?\n\n',['$1' ...
    'pcoupl          = Berendsen' newline ...
    'pcoupltype      = ' Settings.Isotropy newline ...
    'tau-p           = ' num2str(tau_p) newline ...
    'nstpcouple      = ' num2str(nstpcouple) newline ...
    'compressibility = ' Compressibility newline ...
    'ref-p           = ' Target_P newline newline]);

% MDP_Template = regexprep(MDP_Template,'(\; Pressure Coupling).+?\r','$1');
% MDP_Template = regexprep(MDP_Template,'(pcoupl += *)(.+?)( *);','$1Berendsen$3;');
% MDP_Template = regexprep(MDP_Template,'(pcoupltype += *)(.+?)( *);',['$1' Settings.Isotropy '$3;']);
% MDP_Template = regexprep(MDP_Template,'(tau-p += *)(.+?)( *);',['$1 ' num2str(tau_p) '$3;']);
% MDP_Template = regexprep(MDP_Template,'(nstpcouple += *)(.+?)( *);',['$1 ' num2str(nstpcouple) '$3;']);
% MDP_Template = regexprep(MDP_Template,'(compressibility += *)(.+?)( *);',['$1 ' Compressibility '$3;']);
% MDP_Template = regexprep(MDP_Template,'(ref-p += *)(.+?)( *);',['$1 ' Target_P '$3;']);
    
% Pair it with velocity rescale thermostat + small time constant
MDP_Template = regexprep(MDP_Template,'(\; Temperature Coupling\n).+?\n\n',['$1' ...
    'tcoupl     = v-rescale' newline ...
    'tc-grps    = System' newline ...
    'tau-t      = ' num2str(tau_t) newline ...
    'nsttcouple = ' num2str(nsttcouple) newline ...
    'ref-t      = ' Target_T newline newline]);

% MDP_Template = regexprep(MDP_Template,'(tcoupl += *)(.+?)( +);','$1v-rescale$3;');
% MDP_Template = regexprep(MDP_Template,'(tau-t += *)(.+?)( +);',['$1 ' num2str(tau_t) '$3;']);
% MDP_Template = regexprep(MDP_Template,'(nsttcouple += *)(.+?)( +);',['$1 ' num2str(nsttcouple) '$3;']);
    
% Save MDP file, topology file can be reused
MDP_Filename = fullfile(Settings.WorkDir,'Equil_System.mdp');
fidMDP = fopen(MDP_Filename,'wt');
fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
fclose(fidMDP);

TPR_File = fullfile(Settings.WorkDir,'Equil_System.tpr');
MDPout_File = fullfile(Settings.WorkDir,'Equil_System_out.mdp');
GrompLog_File = fullfile(Settings.WorkDir,'Equil_System_Grompplog.log');

% Only applies to polarizable models
ndx_filename = fullfile(Settings.WorkDir,'Equil_System.ndx');
ndx_add = add_polarization_shells(Settings,Settings.SuperCellFile,...
    'ndx_filename',ndx_filename,'add_shells',false);

FEquil_Grompp = [Settings.gmx_loc Settings.grompp ' -c ' windows2unix(Settings.SuperCellFile) ...
    ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Settings.Topology_File) ...
    ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
    ndx_add ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GrompLog_File)];
[state,~] = system(FEquil_Grompp);
    
% Catch errors in grompp
if state ~= 0
    error(['Error running GROMPP. Problem command: ' newline FEquil_Grompp]);
else
    delete(GrompLog_File)
end
    
% Prepare Equilibration mdrun command
Log_File = fullfile(Settings.WorkDir,'Equil_System.log');
Energy_file = fullfile(Settings.WorkDir,'Equil_System.edr');
TRR_File = fullfile(Settings.WorkDir,'Equil_System.trr');
Equilibrated_Geom_File = fullfile(Settings.WorkDir,['Equil_System.' Settings.CoordType]);
    
mdrun_command = [Settings.gmx Settings.mdrun ' -s ' windows2unix(TPR_File) ...   
    ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
    ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Equilibrated_Geom_File) ...
    ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Equil_System')) ...
    Settings.mdrun_opts];

if Settings.Table_Req
    mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
end
    
% Run system Equilibration
if Settings.Verbose
    disp(['Beginning System Pre-Equilibration for ' num2str(Settings.PreEquilibration) ' ps...'] )
end
mintimer = tic;
[state,mdrun_output] = system(mdrun_command);
if state == 0
    if Settings.Verbose
        disp(['System Successfully Equilibrated! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
    end
    while true
        try
            copyfile(Equilibrated_Geom_File,Settings.SuperCellFile)
            break
        catch
            pause(10)
        end
    end
else
    try % Clean up
        [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
        [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.OuterDir) ' -iname "*core*" ' Settings.pipe ' xargs rm']);
    catch me
        if Settings.Verbose
            disp(me.message)
        end
    end
    Settings = Inp_Settings;
	if Settings.MDP.dt/2 >= Settings.MinTimeStep
        if Settings.Verbose
            disp('Equilibration failed. Reducing time step.')
        end
        %Settings.QECompressibility = Settings.QECompressibility_init;
        Settings.MDP.dt = Settings.MDP.dt/2;
        Settings.Output_Coords = Settings.Output_Coords*2;
        save(fullfile(Directory,'TempJobInfo.mat'),'Settings');
        Output = MD_PreEquilibrate(Directory);
        return
	else
        if Settings.Verbose
            disp('Equilibration failed. Reducing time step did not resolve.')
            disp(mdrun_output);
            disp(['Error running mdrun for system pre-equilibration. Problem command: ' newline mdrun_command]);
        end
        Output.Aborted = true;
        return
	end
end

%     system(['wsl source ~/.bashrc; echo "5 8 15 0" ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o ' windows2unix(strrep(Energy_file,'.edr','.xvg'))])
%     system(['echo "5 8 15 0" | gmx_d energy -f ' Energy_file ' -o ' strrep(Energy_file,'.edr','.xvg')])
%     En_xvg_file = fullfile(Settings.WorkDir,'Equil_System.xvg');
%     Data = import_xvg(En_xvg_file);
%     plot(Data(:,1),Data(:,2)./Settings.N_atoms) % Potential
%     plot(Data(:,1),Data(:,3)./Settings.N_atoms) % Conversved Energy
%     plot(Data(:,1),(10^3).*Data(:,4)./Settings.N_atoms) % Volume

% Update MDP file and overwrite
MDP_Template = regexprep(Settings.MDP_Template,'(gen-vel *= *)yes(.+?\n)','$1no $2');
MDP_Template = strtrim(regexprep(MDP_Template,'gen-temp.+?generation',''));
fidMDP = fopen(Settings.MDP_in_File,'wt');
fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
fclose(fidMDP);

% Rerun Grompp
if isfile(Settings.Traj_Conf_File)
    delete(Settings.Traj_Conf_File)
end
GROMPP_command = [Settings.gmx_loc Settings.grompp ' -c ' windows2unix(Settings.SuperCellFile) ...
    ' -f ' windows2unix(Settings.MDP_in_File) ' -p ' windows2unix(Settings.Topology_File) ...
    ' -o ' windows2unix(Settings.Traj_Conf_File) ' -po ' windows2unix(Settings.MDP_out_File) ...
    ndx_add ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(Settings.GrompLog_File)];
[errcode,~] = system(GROMPP_command);

% Catch error in grompp
if errcode ~= 0
    error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
end

if Settings.Verbose
    disp('*** Pre-Equilibration of System Complete ***')
end

% Remove backups
if Settings.Delete_Backups
    system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
        ' -iname "#*#" -delete']);
end

end