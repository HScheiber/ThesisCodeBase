function Completed_Elements = GD_Innerloop(Salt,Structures,datadir,ReCalc_Salts,ReCalc_RefAtom,ReCalc_BSSE,Completed_Elements)
total_timer = tic;
disp(['Begin Data Scraping for: ' Salt])

[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Get Lone ion energies and Basis set superposition error
% Metal
if ReCalc_RefAtom && ~contained_in_cell(Metal,Completed_Elements) && ~isempty(Metal) % Only needs to be done once.
    % Neutral atom energies
    timer = tic;
    disp(['Extracting Atom Energies for: ' Metal '...'])
    Data = Extract_CRYSTAL_Energies(Metal,'Atom');
    filename = [Metal '_Atom_Total_Energies.mat'];
    parsave([datadir filesep filename],Data,'Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Ion energies
    timer = tic;
    disp(['Extracting Ion Energies for: ' Metal '...'])
    Data = Extract_CRYSTAL_Energies(Metal,'Ion');
    filename = [Metal '_Ion_Total_Energies.mat'];
    parsave([datadir filesep filename],Data,'Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
end
if ReCalc_BSSE && ~contained_in_cell(Metal,Completed_Elements) && ~isempty(Metal) % Only needs to be done once.
    % BSSE Energies
    timer = tic;
    disp(['Extracting BSSE Correction Energies for: ' Metal '...'])
    filename = [Metal '_BSSE_Energies.mat'];
    Metal_BSSE_Data = Extract_CRYSTAL_BSSE(Metal,Data);
    parsave([datadir filesep filename],Metal_BSSE_Data,'Metal_BSSE_Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    Completed_Elements{end+1} = Metal;
elseif ~isempty(Metal) % Otherwise load  BSSE Corrections from memory
    timer = tic;
    disp(['Loading BSSE Correction Energies for: ' Metal '...'])
    filename = [Metal '_BSSE_Energies.mat'];
    X = load([datadir filesep filename],'Metal_BSSE_Data');
    Metal_BSSE_Data = X.Metal_BSSE_Data;
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
else
    Metal_BSSE_Data = [];
end

% Halide 
if ReCalc_RefAtom && ~contained_in_cell(Halide,Completed_Elements) && ~isempty(Halide) % Only needs to be done once.
    % Neutral atom energies
    timer = tic;
    disp(['Extracting Atom Energies for: ' Halide '...'])
    Data = Extract_CRYSTAL_Energies(Halide,'Atom');
    filename = [Halide '_Atom_Total_Energies.mat'];
    parsave([datadir filesep filename],Data,'Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Ion energies
    timer = tic;
    disp(['Extracting Ion Energies for: ' Halide '...'])
    Data = Extract_CRYSTAL_Energies(Halide,'Ion');
    filename = [Halide '_Ion_Total_Energies.mat'];
    parsave([datadir filesep filename],Data,'Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
end

if ReCalc_BSSE && ~contained_in_cell(Halide,Completed_Elements) && ~isempty(Halide)
    % BSSE Energies
    timer = tic;
    disp(['Extracting BSSE Correction Energies for: ' Halide '...'])
    filename = [Halide '_BSSE_Energies.mat'];
    Halide_BSSE_Data = Extract_CRYSTAL_BSSE(Halide,Data);
    parsave([datadir filesep filename],Halide_BSSE_Data,'Halide_BSSE_Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    Completed_Elements{end+1} = Halide;
elseif ~isempty(Halide) % Otherwise load  BSSE Corrections from memory
    timer = tic;
    disp(['Loading BSSE Correction Energies for: ' Halide '...'])
    filename = [Halide '_BSSE_Energies.mat'];
    X = load([datadir filesep filename],'Halide_BSSE_Data');
    Halide_BSSE_Data = X.Halide_BSSE_Data;
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
else
    Halide_BSSE_Data = [];
end

%% Scrape the Energies of requested structure types
N = length(Structures);
for i=1:N
	Structure = Structures{i};
    
    % Total energy PES
    timer = tic;
    if ReCalc_Salts
        disp(['Extracting total energies for ' Structure ' ' Salt '...'])
        Labend = Label_replace(Structure);
        Label = [Salt Labend];
        directory = ([Salt filesep Structure]);
        Data = Extract_CRYSTAL_PES(Label,directory,Structure,false);
        if ~isempty(Data)
            filename = [Salt '_' Structure '_Total_Energies.mat'];
            parsave([datadir filesep filename],Data,'Data')
        end
    else
        disp(['Loading total energies for ' Structure ' ' Salt '...'])
        Labend = Label_replace(Structure);
        Label = [Salt Labend];
        filename = [Salt '_' Structure '_Total_Energies.mat'];
        X = load([datadir filesep filename],'Data');
        Data = X.Data;
    end
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Lattice energy PES
    timer = tic;
    disp(['Calculating lattice energies for ' Structure ' ' Salt '...'])
    filename = [Salt '_' Structure '_Lattice_Energies.mat'];
    Data = Generate_Lattice_energy_CRYSTAL(Salt,Label,Structure,datadir,'Lattice');
    if ~isempty(Data)
        parsave([datadir filesep filename],Data,'Data')
    end
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Subtract BSSE
    timer = tic;
    disp(['Calculating BSSE corrected lattice energies for ' Structure ' ' Salt '...'])
    filename = [Salt '_' Structure '_BSSE_Corrected_Lattice_Energies.mat'];
    if ~isempty(Data)
        BSSE_Corrected_Data = Lattice_Energy_BSSE_Correction(Data,...
            Metal_BSSE_Data,Halide_BSSE_Data,Label);
        parsave([datadir filesep filename],BSSE_Corrected_Data,'BSSE_Corrected_Data')
    end
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Cohesive energy PES
    timer = tic;
    disp(['Calculating cohesive energies for ' Structure ' ' Salt '...'])
    filename = [Salt '_' Structure '_Cohesive_Energies.mat'];
    Data = Generate_Lattice_energy_CRYSTAL(Salt,Label,Structure,datadir,'Cohesive');
    if ~isempty(Data)
        parsave([datadir filesep filename],Data,'Data')
    end
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    % Subtract BSSE
    timer = tic;
    disp(['Calculating BSSE corrected cohesive energies for ' Structure ' ' Salt '...'])
    filename = [Salt '_' Structure '_BSSE_Corrected_Cohesive_Energies.mat'];
    if ~isempty(Data)
        BSSE_Corrected_Data = Lattice_Energy_BSSE_Correction(Data,...
            Metal_BSSE_Data,Halide_BSSE_Data,Label);
    end
    parsave([datadir filesep filename],BSSE_Corrected_Data,'BSSE_Corrected_Data')
    elTime = toc(timer);
    disp(['Done (' num2str(elTime) ' s).'])
    
    disp('-------------------------------------------------------------')
    disp(['Data collection complete for ' Salt ' ' Structure ' (' num2str(i) '/' num2str(N) ').'])
    disp('-------------------------------------------------------------')
end

%% Display Total Time
elapsedTime = toc(total_timer);
if elapsedTime < 60
    disp('#############################################################')
    disp(['Data collection complete for ' Salt ' Time Elapsed: ' num2str(elapsedTime) ' s.'])
    disp('#############################################################')
    disp(newline)
else
    elapsedMin = floor(elapsedTime/60);
    elapsedSec = elapsedTime-(elapsedMin*60);
    disp('#############################################################')
    disp(['Data collection complete for ' Salt ' Time Elapsed: ' num2str(elapsedMin) ' m ' num2str(elapsedSec) ' s.'])
    disp('#############################################################')
    disp(newline)
end
end