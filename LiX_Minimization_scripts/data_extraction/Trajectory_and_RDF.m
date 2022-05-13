% Trajectory_and_RDF
% Data analysis script for GROMACS MD trajectories. Makes RDFs and compact
% trajectories if requested.
function Data = Trajectory_and_RDF(Directory,~)

%% OPTIONS
Generate_Compressed_Trj = true; % If true, generates a compressed .xtc trajectory file
CompTrj_pbc = 'atom'; % Option -pbc sets the type of periodic boundary condition treatment
CompTrj_unitcell = 'tric'; % One of 'rect', 'compact', 'tric'. Sets the unit cell representation
CompTrj_dt = 10; % [ps] Change time step to larger one. 0 for default

Generate_RDF = true; % If true, generate RDF data
RDF_Breakpoints = 0:1000:10000;
RDF_dt = 100;% [ps] Change time step to larger one. 0 for default

% End of options
old_directory = pwd;
cd(Directory)

Data = struct;

if ispc % for testing
    gmx = 'wsl source ~/.bashrc; gmx_d';
    sys = @(inp) system(inp);
elseif isunix
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
    if strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'sea')
        gmx = 'gmx_d';
        sys = @(inp) system(inp);
    elseif strcmpi(Server,'bel') || strcmpi(Server,'blg')
        gmx = 'gmx_d';
        sys = @(inp) system_def(inp); % Needed to circumvent error
    elseif strcmpi(Server,'pat')
        gmx = 'gmx_d';
        sys = @(inp) system(inp); 
    end
else
    error('Unknown machine.')
end

% Find the trr (trajectory) file
Files = dir('*.trr');

% Get filebase from trr name
[~,Filebase,~] = fileparts(Files.name);

Salt = regexp(Filebase,'.+?(?=_)','match','ONCE');
[Metal,Halide] = Separate_Metal_Halide(Salt);

% Assign files
TRR_File = [Filebase '.trr']; % trajectory output file
TPR_File = [Filebase '.tpr']; % Trajectory input file
XTC_File = [Filebase '.xtc']; % Compressed Trajectory output file
EDR_File = [Filebase '.edr']; % Binary energy file
Structure_File = [Filebase '.g96']; % XVG energy file for data analysis
MAT_File = [Filebase '.mat']; % Contains info about number of particles
XVG_Energy_File = [Filebase '_Energies.xvg']; % XVG energy file for data analysis

if Generate_Compressed_Trj && ~exist(XTC_File,'file')
    Compress_Trj_Cmd = ['echo 0 | ' gmx ' trjconv -f ' TRR_File ' -s ' TPR_File ...
        ' -pbc ' CompTrj_pbc ' -ur ' CompTrj_unitcell ' -o ' XTC_File ' -dt ' num2str(CompTrj_dt)];
    sys(Compress_Trj_Cmd);
end

%% Generate total energy vs time, total V vs time, pressure vs time, T vs time.
Energy_Command = ['echo 4 8 9 13 0 | ' gmx ' energy -f ' EDR_File ' -o ' XVG_Energy_File];
sys(Energy_Command);

Data.Energy = import_xvg(XVG_Energy_File);

MatData = load(MAT_File);
Data.Energy(:,2) = Data.Energy(:,2)*2./MatData.N_total;
Data.Energy(:,5) = Data.Energy(:,5)*2./MatData.N_total;

%% Generate RDF
if Generate_RDF
    N = length(RDF_Breakpoints);

    rdf_tic = tic;
    disp(['Generating Initial RDFs (1/' num2str(N) ')']);
    
    MM_RDF_Command = [gmx ' rdf -f ' Structure_File ' -s ' TPR_File ' -o RDF.xvg' ...
        ' -ref ' Metal ' -sel ' Metal ' -dt ' num2str(RDF_dt)];
    sys(MM_RDF_Command);
    Data.RDF.MM = cell(N,2);
    Data.RDF.MM{1,1} = 0;
    Data.RDF.MM{1,2} = import_xvg('RDF.xvg');
    delete RDF.xvg

    HH_RDF_Command = [gmx ' rdf -f ' Structure_File ' -s ' TPR_File ' -o RDF.xvg' ...
        ' -ref ' Halide ' -sel ' Halide ' -dt ' num2str(RDF_dt)];
    sys(HH_RDF_Command);
    Data.RDF.HH = cell(N,2);
    Data.RDF.HH{1,1} = 0;
    Data.RDF.HH{1,2} = import_xvg('RDF.xvg');
    delete RDF.xvg

    MH_RDF_Command = [gmx ' rdf -f ' Structure_File ' -s ' TPR_File ' -o RDF.xvg' ...
        ' -ref ' Metal ' -sel ' Halide ' -dt ' num2str(RDF_dt)];
    sys(MH_RDF_Command);
    Data.RDF.MH = cell(N,2);
    Data.RDF.MH{1,1} = 0;
    Data.RDF.MH{1,2} = import_xvg('RDF.xvg');
    delete RDF.xvg

    rdf_toc = duration(seconds(toc(rdf_tic)),'Format','hh:mm:ss');     
    disp(['Initial RDFs Complete (1/' num2str(N) '). Time Elapsed: ' char(rdf_toc)]);
    for i=2:N
        rdf_tic = tic;
        Initial_Pnt = num2str(RDF_Breakpoints(i-1));
        Final_Pnt = num2str(RDF_Breakpoints(i));
        disp(['Generating RDFs from ' Initial_Pnt ' - ' Final_Pnt ' ps. (' num2str(i) '/' num2str(N) ')']);
        

        MM_RDF_Command = [gmx ' rdf -f ' XTC_File ' -s ' TPR_File ' -o RDF.xvg' ...
            ' -ref ' Metal ' -sel ' Metal ' -b ' Initial_Pnt ' -e ' Final_Pnt ' -dt ' num2str(RDF_dt)];
        [ErrorCode,~] = sys(MM_RDF_Command);
        if ErrorCode ~= 0
            rdf_toc = duration(seconds(toc(rdf_tic)),'Format','hh:mm:ss');   
            disp(['RDFs ' Initial_Pnt ' - ' Final_Pnt ' ps failed. Qutting. Time Elapsed: ' char(rdf_toc)])
            break
        end
        Data.RDF.MM{i,1} = RDF_Breakpoints(i);
        Data.RDF.MM{i,2} = import_xvg('RDF.xvg');
        delete RDF.xvg

        HH_RDF_Command = [gmx ' rdf -f ' XTC_File ' -s ' TPR_File ' -o RDF.xvg' ...
            ' -ref ' Halide ' -sel ' Halide ' -b ' Initial_Pnt ' -e ' Final_Pnt ' -dt ' num2str(RDF_dt)];
        [ErrorCode,~] = sys(HH_RDF_Command);
        if ErrorCode ~= 0
            rdf_toc = duration(seconds(toc(rdf_tic)),'Format','hh:mm:ss');   
            disp(['RDFs ' Initial_Pnt ' - ' Final_Pnt ' ps failed. Qutting. Time Elapsed: ' char(rdf_toc)])
            break
        end
        Data.RDF.HH{i,1} = RDF_Breakpoints(i);
        Data.RDF.HH{i,2} = import_xvg('RDF.xvg');
        delete RDF.xvg

        MH_RDF_Command = [gmx ' rdf -f ' XTC_File ' -s ' TPR_File ' -o RDF.xvg' ...
            ' -ref ' Metal ' -sel ' Halide ' -b ' Initial_Pnt ' -e ' Final_Pnt ' -dt ' num2str(RDF_dt)];
        [ErrorCode,~] = sys(MH_RDF_Command);
        if ErrorCode ~= 0
            rdf_toc = duration(seconds(toc(rdf_tic)),'Format','hh:mm:ss');   
            disp(['RDFs ' Initial_Pnt ' - ' Final_Pnt ' ps failed. Qutting. Time Elapsed: ' char(rdf_toc)])
            break
        end
        Data.RDF.MH{i,1} = RDF_Breakpoints(i);
        Data.RDF.MH{i,2} = import_xvg('RDF.xvg');
        delete RDF.xvg
        
        rdf_toc = duration(seconds(toc(rdf_tic)),'Format','hh:mm:ss');   
        disp(['RDFs ' Initial_Pnt ' - ' Final_Pnt ' ps complete (' ...
            num2str(i) '/' num2str(N) '). Time Elapsed: ' char(rdf_toc)]);
    end
end

% Return
cd(old_directory);

end
