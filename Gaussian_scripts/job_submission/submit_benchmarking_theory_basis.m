function submit_benchmarking_theory_basis(Label)
%% Inputs
hours = 3; % Max time for job (hours)
mins = 0; % Max time for job (minutes)

pausetime = 1; % time to pause
nMols_per_Task = -1; % -1 to fix the number of cores
nCores = 32;
nTasks_per_Node = 32;
mempernode = "-1"; % Memory request for server (default = "-1")
dispersion = ""; % Set to "GD3BJ" / "GD2" / "PFD" / "GD3" / "" (no dispersion);
conv = 8; % Convergence criteria
memory = "64000MB";
guess = ""; % Leave as "" for default

% Theory
% theory = {"B3LYP","B2PLYP","B2PLYPD3","LC-LSDA","LC-PBEPBE","wB97XD","DSDPBEP86",...
%     "HSEH1PBE","wB97X","CAM-B3LYP","LC-wHPBE","B2GPPLYP","PW91PW91","mPW1PW91",...
%     "PW1PW","MP4","CISD","CCSD(T)"};
% theory = {"LC-LSDA","LSDA","PBEPBE","PBE1PBE","DSDPBEP86","B2GPPLYP","HSEH1PBE","wB97X","wB97","B3LYP",...
%     "CAM-B3LYP","LC-wHPBE","B2PLYP","PW91PW91","mPW1PW91","PW1PW","TPSSTPSS","RevTPSSRevTPSS","HISSbPBE"}; 
 theory = {"LSDA", "PBEPBE", "RevTPSSRevTPSS"};

% Basis set
%basis = {"6-311+G*","CBSB7++","AUG-cc-pVDZ","AUG-cc-pVTZ","AUG-cc-pVQZ","Def2SVPP","Def2TZVPP","Def2QZVPP"};
%basis = {"pob-TZVP","7-311G-jansen","7-311G-crystal","7-311G-nada"};
basis = {"Def2TZVP"};
%% Code
% Autochoose home folder based on system type
if ispc
    home = "C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase"; % PC
elseif isunix
    home = "/home/scheiber/ThesisCodeBase"; % Cedar/Graham
else
    home = input("Please input home bin directory.\n","s");
end

% Autochoose template file
if Label == "F"
    template_file = "F_Ion.template";
elseif Label == "FF"
    template_file = "FF_Ion_Pair.template";
elseif Label == "LiF"
    template_file = "LiF_Ion_Pair.template";
elseif Label == "LiLi"
    template_file = "LiLi_Ion_Pair.template";
elseif Label == "Li"
    template_file = "Li_Ion.template";
elseif Label == "Cl"
    template_file = "Cl_Ion.template";
elseif Label == "Br"
    template_file = "Br_Ion.template";
elseif Label == "I"
    template_file = "I_Ion.template";
else
    error("Unknown label: """ + Label +...
		""". Choose one of ""F"", ""FF"", ""LiF"", ""LiLi"", ""Li"", ""Cl"", ""Br"", ""I"".")
end

% For effective core potentials
if Label == "I"
    gen = "GenECP";
else
    gen = "Gen";
end

nTasks = num2str(nCores);
hours_calc = num2str(hours);
mins_calc = num2str(mins);
nCores_per_Node = num2str(nTasks_per_Node);
ECPlink = "";
workDir = string(pwd);
submit_Gaussian = "submit_Gaussian.pl";

for i = 1:length(theory)
    
    l1s = theory{i};  
    l1d = Label + "_" + regexprep(theory{i},"\((.?)\)","_$1");
    
    % Set additional environment variables option for submission script
    if (l1s == "DSDPBEP86") && (dispersion == "GD3BJ")
        ExtraEV = "0";
    elseif (l1s == "B2GPPLYP") && (dispersion == "GD3BJ")
        ExtraEV = "1";
    else
        ExtraEV = "-1";
    end
    
    % create/change dir
    if ~exist(l1d,"dir")
        mkdir(l1d);
    end
    
    dir1 = workDir + filesep + l1d;
    cd(dir1)

    for j = 1:length(basis)

        l2s = basis{j};
        l2d = "BASS_" + basis{j};

        % create/change dir
        if ~exist(l2d,"dir")
            mkdir(l2d);
        end
        dir2 = workDir + filesep + l1d + filesep + l2d;
        cd(dir2)

        % Open template
        fid = fopen( home + filesep + "templates" + filesep + template_file,"rt");
        X = fread(fid);
        fclose(fid);
        newtext = char(X.');

        % Replace strings
        % theory
        if l1s == "B2GPPLYP"
            newtext = strrep(newtext,"##THERT##","B2GP-PLYP");
            newtext = strrep(newtext,"##THER##","b2plyp");
            newtext = strrep(newtext,"##IOPS##",...
                newline + "iop(3/125=0360003600,3/78=0640006400,3/76=0350006500,3/77=1000010000)");
            if (dispersion ~= "GD3BJ") && (dispersion ~= "")
                error("Unknown parameters for dispersion type " +...
                    dispersion + " for " + l1s);
            end
        elseif l1s == "DSDPBEP86"
            newtext = strrep(newtext,"##THERT##","DSD-PBEP86");
            newtext = strrep(newtext,"##THER##","b2plyp");
            newtext = strrep(newtext,"##IOPS##",...
                newline + "iop(3/125=0220005200,3/78=0440004400,3/76=0310006900,3/74=1004)");
            if (dispersion ~= "GD3BJ") && (dispersion ~= "")
                error("Unknown parameters for dispersion type " +...
                    dispersion + " for " + l1s);
            end
        elseif l1s == "PW1PW"
            newtext = strrep(newtext,"##THERT##","PW1PW");
            newtext = strrep(newtext,"##THER##","PW91PW91");
            newtext = strrep(newtext,"##IOPS##",...
                newline + "iop(3/76=1000002000,3/77=0800000000,3/78=1000000000)");
            if (dispersion ~= "")
                error("Unknown parameters for dispersion type " +...
                    dispersion + " for " + l1s);
            end
        else
            newtext = strrep(newtext,"##THERT##",l1s);
            newtext = strrep(newtext,"##THER##",l1s);
            newtext = strrep(newtext,"##IOPS##","");
        end
        
        % Basis set
        if l2s == "pob-TZVP"
            newtext = strrep(newtext,"##BASS##","Gen");
            newtext = strrep(newtext,"##BASST##","pob-TZVP");
            % Custom basis sets
            newtext = strrep(newtext,"##GENBAS##","@" + home + filesep +...
                "basis_sets" + filesep + "pob_TZVP.gbs" + newline + newline);
            if Label == "I"
                ECPlink = "@" + home + filesep + "basis_sets" + filesep +...
                    "I_ECP_pob_TZVP.ecp" + newline + newline;
            end
        elseif l2s == "7-311G-crystal"
            newtext = strrep(newtext,"##BASS##","Gen");
            newtext = strrep(newtext,"##BASST##","7-311G-crystal");
            % Custom basis sets
            newtext = strrep(newtext,"##GENBAS##","@" + home + filesep +...
                "basis_sets" + filesep + "7_311G_crystal.gbs" + newline + newline);
        elseif l2s == "7-311G-jansen"
            newtext = strrep(newtext,"##BASS##","Gen");
            newtext = strrep(newtext,"##BASST##","7-311G-jansen");
            % Custom basis sets
            newtext = strrep(newtext,"##GENBAS##","@" + home + filesep +...
                "basis_sets" + filesep + "7_311G_jansen.gbs" + newline + newline);
        elseif l2s == "7-311G-nada"
            newtext = strrep(newtext,"##BASS##","Gen");
            newtext = strrep(newtext,"##BASST##","7-311G-nada");
            % Custom basis sets
            newtext = strrep(newtext,"##GENBAS##","@" + home + filesep +...
                "basis_sets" + filesep + "7_311G_nada.gbs" + newline + newline);
        else
            newtext = strrep(newtext,"##BASS##", l2s);
            newtext = strrep(newtext,"##BASST##", l2s);
            newtext = strrep(newtext,"##GENBAS##","");
        end
        
        newtext = strrep(newtext,"##CORES##",nTasks);
        newtext = strrep(newtext,"##CONV##", num2str(conv));
        newtext = strrep(newtext,"##MEMORY##", num2str(memory));
        newtext = strrep(newtext,newline + "##GENECP##",ECPlink);
        
        % Dispersion correction
        if dispersion == ""
            newtext = strrep(newtext,"##DISPT##","");
            newtext = strrep(newtext,"##DISP##","");
        else
            newtext = strrep(newtext,"##DISP##","EmpiricalDispersion=" +...
                dispersion);
            newtext = strrep(newtext,"##DISPT##",...
                "and " + dispersion + " empirical dispersion correction");
        end
		
		% Modify guess type
		if guess == ""
			newtext = strrep(newtext,"##GUESS##","");
		else
			newtext = strrep(newtext,"##GUESS##"," Guess=" + guess);
		end
        
        % Save file
        fid2 = fopen( l1d + "_" + l2s + ".com","wt");
        fwrite(fid2,newtext);
        fclose(fid2);

        % Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
        command = submit_Gaussian + " -1 " + nTasks + " " + nCores_per_Node +...
        " -1 " + mempernode + " " + hours_calc + " " + mins_calc + " " + l1d + "_" + l2s +...
        " " + l1d + "_" + l2s + ".com " + l1d + "_" + l2s + " 1 " + ExtraEV;
        disp(l1d + " " + l2s)

        if ~ispc
            system(command);
            pause(pausetime);
        end

        cd(dir1)
    end
    cd(workDir)
end
end