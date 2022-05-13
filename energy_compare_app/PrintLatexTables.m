function LoadbarObj = PrintLatexTables(LoadbarObj,Struct_Data,Ymin,filename,Data_Types,Basis_Set)

    Remove_D3gCP = false; % Set true if you wish to remove D3 and gCP markers in the output
    Structure_List = {'Rocksalt' 'Wurtzite' 'BetaBeO' 'CsCl' 'FiveFive' 'NiAs' 'Sphalerite'};
    N_SL = length(Structure_List);

    N_Col = size(Struct_Data,1);
    
    Species = Struct_Data(:,4);
    % Rename species to order correctly
    Species = strrep(Species,'LiF','AAA');
    Species = strrep(Species,'LiCl','BBB');
    Species = strrep(Species,'LiBr','CCC');
    Species = strrep(Species,'LiI','DDD');
    Species = strrep(Species,'NaCl','EEE');
    
    Structure_types = regexprep(Struct_Data(:,2),'(Rocksalt|Wurtzite)','AAA$1');
    Theory_types = regexprep(Struct_Data(:,1),'_','-');
    Theory_types = regexprep(Theory_types,'(JC|TF)','ZZZ$1'); %% Make sure empirical are sorted last
    Theory_types = strrep(Theory_types,'Experiment','AAA'); %% Make sure experiment are sorted first
    Energy_min = cell2mat(Ymin(:));
    a_min = cell2mat(Struct_Data(:,3));
    b_min = cell2mat(Struct_Data(:,5));
    c_min = cell2mat(Struct_Data(:,6));
    Volume_min = cell2mat(Struct_Data(:,7));
    BL_min = cell2mat(Struct_Data(:,8));
    alpha_min = cell2mat(Struct_Data(:,10));
    beta_min = cell2mat(Struct_Data(:,11));
    gamma_min = cell2mat(Struct_Data(:,12));
    aT_min = cell2mat(Struct_Data(:,14));
    bT_min = cell2mat(Struct_Data(:,15));
    cT_min = cell2mat(Struct_Data(:,16));
    
    FC_min = cell(N_Col,1);
    Opt_Type = cell(N_Col,1);
    
    if length(Data_Types) > 1
        warning('Multiple data types chosen. Assuming fully optimized.')
        Sectionheader = 'for fully optimized unit cells with fixed symmetry.';
        Opt_Type = 'from Full Geometry Optimization.';
        Label = 'Full';
    else
        switch Data_Types{1}
            case 'Rigid Structure'
                Sectionheader = 'for volume-optimized unit cells with fixed fractional coordinates.';
                Opt_Type = 'From Energy Curve.';
                Label = 'Curve';
                Note = '';
            case 'Cell Opt'
                Sectionheader = 'for volume-optimized unit cells with fixed fractional coordinates.';
                Opt_Type = 'with Fixed Fractional Coordinates.';
                Label = 'Volume';
                Note = '$^{\ddagger}$';
            case 'Full Opt'
                Sectionheader = 'for fully optimized unit cells with fixed symmetry.';
                Opt_Type = 'from Full Geometry Optimization.';
                Label = 'Full';
                Note = '';
            case 'Atom Optimized'
                Sectionheader = 'for atom-optimized unit cells with fixed lattice parameters.';
                Opt_Type = 'from Fractional Coordinate Optimization.';
                Label = 'Atom';
                Note = '';
            case 'Full Opt (No Sym)'
                Sectionheader = 'for fully optimized unit cells without any constraints.';
                Opt_Type = 'from Full Geometry Optimization without symmetry.';
                Label = 'Full_SG1';
                Note = '';
            case 'Full Opt (P = 1 atm)'
                Sectionheader = 'for fully optimized unit cells with fixed symmetry and 1 bar external pressure.';
                Opt_Type = 'from Full Geometry Optimization at 1 bar pressure.';
                Label = 'Full_P1';
                Note = '';
        end
    end
    
    for i=1:N_Col
        FC = Struct_Data{i,9}';
        for x=1:size(FC,2)
            FC(2:4,x) = num2cell(mod(round([FC{2:4,x}].*1e4)/1e4,1)');
        end
        
        % Sort by metal first then halide
        FC(1,:) = regexprep(FC(1,:),'L(I|i)','AAA');
        FC(1,:) = regexprep(FC(1,:),'N(A|a)','AAB');
        FC(1,:) = regexprep(FC(1,:),'K','AAC');
        FC(1,:) = regexprep(FC(1,:),'R(B|b)','AAD');
        FC(1,:) = regexprep(FC(1,:),'C(S|s)','AAE');
        [~,SortIdx] = sort(FC(1,:));
        FC = FC(:,SortIdx);
        FC(1,:) = strrep(FC(1,:),'AAA','LI');
        FC(1,:) = strrep(FC(1,:),'AAB','Na');
        FC(1,:) = strrep(FC(1,:),'AAC','K');
        FC(1,:) = strrep(FC(1,:),'AAD','RB');
        FC(1,:) = strrep(FC(1,:),'AAE','CS');
        
        FC = sprintf(' (%5.4f, %5.4f, %5.4f)                           &',FC{2:end,:});
        if ~isempty(Note)
            FC = strrep(FC,')             ',[')' Note]);
        end
        FC_min{i} = FC(1:end-1);
    end
    
    waitbar(0.25,LoadbarObj,'Generating Table...')

    % Merge into table
    Data_Table = table(Species,Structure_types,Theory_types,...
        Energy_min,a_min,b_min,c_min,alpha_min,beta_min,gamma_min,...
        Volume_min,BL_min,FC_min,aT_min,bT_min,cT_min,...
        'VariableNames',{'Species' 'Structure' 'Theory' ...
        'E_min_kJmol' 'a_Ang' 'b_Ang' 'c_Ang' 'alpha_Degree' 'beta_Degree' 'gamma_Degree' ...
        'V_Ang3' 'BL_Ang' 'Frac_Coords' 'aT_Ang' 'bT_Ang' 'cT_Ang'});
    
    % Sort Table
    Data_Table = sortrows(Data_Table,{'Species' 'Theory' 'Structure'},{'ascend','ascend','ascend'});   
    Data_Table2 = sortrows(Data_Table,{'Species' 'Structure' 'Theory'},{'ascend','ascend','ascend'});
    
    % Rename stuff to correct names
    Data_Table{:,1} = strrep(Data_Table{:,1},'AAA','LiF');
    Data_Table{:,1} = strrep(Data_Table{:,1},'BBB','LiCl');
    Data_Table{:,1} = strrep(Data_Table{:,1},'CCC','LiBr');
    Data_Table{:,1} = strrep(Data_Table{:,1},'DDD','LiI');
    Data_Table{:,1} = strrep(Data_Table{:,1},'EEE','NaCl');
    Data_Table{:,2} = strrep(Data_Table{:,2},'AAA','');
    Data_Table{:,3} = strrep(Data_Table{:,3},'AAA','Experiment');
    Data_Table{:,3} = strrep(Data_Table{:,3},'ZZZ','');
    
    % Rename stuff to correct names
    Data_Table2{:,1} = strrep(Data_Table2{:,1},'AAA','LiF');
    Data_Table2{:,1} = strrep(Data_Table2{:,1},'BBB','LiCl');
    Data_Table2{:,1} = strrep(Data_Table2{:,1},'CCC','LiBr');
    Data_Table2{:,1} = strrep(Data_Table2{:,1},'DDD','LiI');
    Data_Table2{:,1} = strrep(Data_Table2{:,1},'EEE','NaCl');
    Data_Table2{:,2} = strrep(Data_Table2{:,2},'AAA','');
    Data_Table2{:,3} = strrep(Data_Table2{:,3},'AAA','Experiment');
    Data_Table2{:,3} = strrep(Data_Table2{:,3},'ZZZ','');
    
    % Find LiI data, specify basis set
    if length(Basis_Set) > 1
        warning('More than one basis set selected. Marking LiI as first basis set in series.')
        Basis_Set = Basis_Set{1};
    end
    
    if strcmp(Basis_Set,'pob-TZVP')
        BS = '$^{*}$';
    elseif strcmp(Basis_Set,'Mpob-TZVP_HECP28MDF-VTZ')
        BS = '$^{\dagger}$';
    else
        BS = '';
    end
    
    for i = 1:size(Data_Table,1)
        if strcmp(Data_Table{i,1}{1},'LiI') && ~contained_in_cell(Data_Table{i,3}{1},{'Experiment' 'JC' 'JC3P' 'JC4P' 'TF'})
            Data_Table{i,3}{1} = [Data_Table{i,3}{1} BS];
        end
    end
    for i = 1:size(Data_Table2,1)
        if strcmp(Data_Table2{i,1}{1},'LiI') && ~contained_in_cell(Data_Table2{i,3}{1},{'Experiment' 'JC' 'JC3P' 'JC4P' 'TF'})
            Data_Table2{i,3}{1} = [Data_Table2{i,3}{1} BS];
        end
    end
    
    %% Energy Tables
    X = load('C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\CRYSTAL\Experimental_Lattice_Energies.mat','Data2');
    ExpData = X.Data2;
    
    % Find the unique list of species
    Q = size(unique(Data_Table(:,2)),1);
    TableAdj = ['{l|l' repmat('l', 1,Q) '}'];
    
    N = size(Data_Table,1);
    
    % Initialize previous salt empty string
    Previous_Salt = '';
    Energy_Table = {};
    
    % Loop through all data entries
    waitbar(0.5,LoadbarObj,'Generating Energy Tables...');
    for Idx = 1:N
        Current_Salt = Data_Table{Idx,1}{1};
        Current_Structure = Data_Table{Idx,2}{1};
        Current_Theory = Data_Table{Idx,3}{1};
        Current_E = num2str(Data_Table{Idx,4},'%6.2f');
        Current_E_num = Data_Table{Idx,4};
        
        [~,Str_Index] = ismember(Current_Structure,Structure_List);
        if Str_Index == 0 
            continue
        end
        
        if strcmp(Current_Structure,'Wurtzite') && strcmp(Current_Theory,'Experiment')
            Current_E = '';
            Current_E_num = NaN;
        elseif strcmp(Current_Structure,'Rocksalt') && strcmp(Current_Theory,'Experiment')
            Err = round(ExpData.(Current_Salt).dE,1,'significant');
            n = 0;
            while (floor(Err*10^n)~=Err*10^n)
                n = n+1;
            end
            Errtxt = num2str(Err,'%1.0e');
            Current_E = [num2str(Current_E_num,['%5.' num2str(n) 'f']) '(' Errtxt(1) ')'];
        end
        
        % Check if salt has changed
        if ~strcmp(Current_Salt,Previous_Salt)
            
            % Previous Energy Table Complete
            if Idx ~= 1
                
                % Finish up last row
                [~,MinIdx] = min(LineEnergies);
                LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];
                
                Favoured = Structure_List{MinIdx};
                Favoured = strrep(Favoured,'Rocksalt','Rock.');
                Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
                Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
                Favoured = strrep(Favoured,'FiveFive','5-5');
                Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
                LinesText{x} = [LinesText{x} '& ' Favoured];
                
                for i = 1:length(LinesData)
                    LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
                    LinesText{x} = [LinesText{x} '&' LinesData{i}];
                end
                
                % Finish title
                for i = 1:N_SL
                    LineTitle = [LineTitle pad(['& \textbf{' Shorten_Structure_Name(Structure_List{i}) '} '],23,'right')];
                end
                
                % Finish up all rows
                LengthT = length(LineTitle);
                LineTitle = [LineTitle '\\ \hline'];
                for i = 1:length(LinesText)
                    if i < length(LinesText)
                        LinesText{i} = ['            ' LinesText{i} '\\'];
                    else
                        LinesText{i} = ['            ' LinesText{i}];
                    end
                end
                
                % Build table
                Energy_Table{end+1} = [Header1 newline ...
                    Header2 newline ...
                    Header3 newline newline ...
                    Line1 newline ...
                    Line2 newline ...
                    Line2b newline ...
                    Line3 newline ...
                    Line4 newline ...
                    LineTitle];
                
                for i = 1:length(LinesText)
                    Energy_Table{end} = [Energy_Table{end} newline ...
                        LinesText{i}];
                end
                
                Energy_Table{end} = [Energy_Table{end} newline ...
                    Line5 newline ...
                    Line6 newline ...
                    Line7 newline ...
                    Line8];
                
            end
            
            x = 0; % keep track of theory number

            % Section headings
            Header1 = ['\subsubsection{' Current_Salt '}'];
            Header2 = ['\label{Optimization_' Label '_' Current_Salt '}'];
            Header3 = ['In this subsection, calculated ' Current_Salt ' properties are reported ' Sectionheader];

            % Start new table
            Line1 = ['%Table: ' Current_Salt ' Energy Results ' Opt_Type];
            Line2 = '\begin{table}[H]';
            Line2b = '    \centering';
            Line3 = '    \begin{adjustwidth}{-0.10in}{-0.10in}% adjust the L and R margins by 0.1 inch';
            Line4 = ['        \begin{tabular}' TableAdj];
            LineTitle = [pad('',46) '& \textbf{Fav.} '];
            Line5 = '        \end{tabular}';
            Line6 = '    \end{adjustwidth}';
            Line7 = ['    \caption{\label{tab:' Current_Salt '_Energy_' Label '} Lattice energies of optimized ' Current_Salt ' crystals with different calculation and crystal structure types. The lowest lattice energy for each calculation is emboldened.}'];
            Line8 = '\end{table}';
            
            Theory_List = {};
            LinesText = {};
        end
        
        if ~contained_in_cell(Current_Theory,Theory_List) % New theory
            
            if strcmp(Current_Salt,Previous_Salt) % Same salt
                [~,MinIdx] = min(LineEnergies);
                LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];
                
                Favoured = Structure_List{MinIdx};
                Favoured = strrep(Favoured,'Rocksalt','Rock.');
                Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
                Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
                Favoured = strrep(Favoured,'FiveFive','5-5');
                Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
                LinesText{x} = [LinesText{x} '& ' Favoured];
                
                for i = 1:length(LinesData)
                    LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
                    LinesText{x} = [LinesText{x} '&' LinesData{i}];
                end
            end
            
            x = x + 1;
            Theory_List{end+1} = Current_Theory;
            if Remove_D3gCP
                Theory_Out = regexprep(Current_Theory,'-D3(-gCP){0,1}','');
            else
                Theory_Out = regexprep(Current_Theory,'(-D3)(-gCP){0,1}','-D3(BJ)$2');
            end
            Theory_Out = regexprep(Theory_Out,'JC(?!.)','JC (SPC/E)');
            Theory_Out = strrep(Theory_Out,'JC3P','JC (TIP3P)');
            Theory_Out = strrep(Theory_Out,'JC4P','JC (TIP4P$_{\text{EW}}$)');
            
            LinesText{x} = pad(['\textbf{' Theory_Out '}'],34,'right');
            LinesData = cell(1,N_SL);
            LineEnergies = nan(1,N_SL);
        end

        LinesData{Str_Index} = Current_E;
        LineEnergies(Str_Index) = Current_E_num;
        Previous_Salt = Current_Salt;
    end
    
    % Build last energy table
    % Finish up last row
    [~,MinIdx] = min(LineEnergies);
    LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];

    Favoured = Structure_List{MinIdx};
    Favoured = strrep(Favoured,'Rocksalt','Rock.');
    Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
    Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
    Favoured = strrep(Favoured,'FiveFive','5-5');
    Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
    LinesText{x} = [LinesText{x} '& ' Favoured];

    for i = 1:length(LinesData)
        LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
        LinesText{x} = [LinesText{x} '&' LinesData{i}];
    end

    % Finish title
    for i = 1:N_SL
        LineTitle = [LineTitle pad(['& \textbf{' Shorten_Structure_Name(Structure_List{i}) '} '],23,'right')];
    end

    % Finish up all rows
    LengthT = length(LineTitle);
    LineTitle = [LineTitle '\\ \hline'];
    for i = 1:length(LinesText)
        if i < length(LinesText)
            LinesText{i} = ['            ' LinesText{i} '\\'];
        else
            LinesText{i} = ['            ' LinesText{i}];
        end
    end
    
    
    % Build table
    Energy_Table{end+1} = [Header1 newline ...
        Header2 newline ...
        Header3 newline newline ...
        Line1 newline ...
        Line2 newline ...
        Line2b newline ...
        Line3 newline ...
        Line4 newline ...
        LineTitle];

    for i = 1:length(LinesText)
        Energy_Table{end} = [Energy_Table{end} newline ...
            LinesText{i}];
    end

    Energy_Table{end} = [Energy_Table{end} newline ...
        Line5 newline ...
        Line6 newline ...
        Line7 newline ...
        Line8];
                
    
    %% Structure Parameter Tables
	waitbar(0.75,LoadbarObj,'Generating Parameter Tables...');
    % Initialize previous salt empty string
    Previous_Salt = '';
    Previous_Structure = '';
    Parameter_Table = {};
    
    y = 0; % Salt counter
    z = 0; % Structure counter
    % Loop through all data entries
    for Idx = 1:N
        Current_Salt = Data_Table2{Idx,1}{1};
        Current_Structure = Data_Table2{Idx,2}{1};
        Current_Theory = Data_Table2{Idx,3}{1};
        Current_a = strrep(num2str(Data_Table2{Idx,5},'%6.4f'),'NaN','');
        Current_b = strrep(num2str(Data_Table2{Idx,6},'%6.4f'),'NaN','');
        Current_c = strrep(num2str(Data_Table2{Idx,7},'%6.4f'),'NaN','');
        
        if ~isnan(Data_Table2{Idx,5}) && ~isnan(Data_Table2{Idx,14})
            Current_a = [Current_a ' (' num2str(Data_Table2{Idx,14},'%6.4f') ')'];
            Current_b = [Current_b ' (' num2str(Data_Table2{Idx,15},'%6.4f') ')'];
            Current_c = [Current_c ' (' num2str(Data_Table2{Idx,16},'%6.4f') ')'];
        end
        Current_alpha = num2str(Data_Table2{Idx,8},'%6.4f');
        Current_beta = num2str(Data_Table2{Idx,9},'%6.4f');
        Current_gamma = num2str(Data_Table2{Idx,10},'%6.4f');
        Current_FC = Data_Table2{Idx,13}{1};
        [Metal,Halide] = Separate_Metal_Halide(Current_Salt);
        Metal = [Metal '$^{+}$'];
        Halide = [Halide '$^{-}$'];
        
        if strcmp(Current_Structure,'Wurtzite') && strcmp(Current_Theory,'Experiment')
            Err_a = round(ExpData.(Current_Salt).Wurtzite.da,1,'significant');
            Err_c = round(ExpData.(Current_Salt).Wurtzite.dc,1,'significant');
            if isnan(Err_a)
                Current_a = '';
            else
                n_a = 0;
                while (floor(Err_a*10^n_a)~=Err_a*10^n_a)
                    n_a = n_a+1;
                end
                Errtxt_a = num2str(Err_a,'%1.0e');
                Current_a = [num2str(ExpData.(Current_Salt).Wurtzite.a,['%5.' num2str(n_a) 'f']) '(' Errtxt_a(1) ')'];
            end
            
            if isnan(Err_c)
                Current_c = '';
            else
                n_c = 0;
                while (floor(Err_c*10^n_c)~=Err_c*10^n_c)
                    n_c = n_c+1;
                end
                Errtxt_c = num2str(Err_c,'%1.0e');
                Current_c = [num2str(ExpData.(Current_Salt).Wurtzite.c,['%5.' num2str(n_c) 'f']) '(' Errtxt_c(1) ')'];
            end
            switch Current_Salt
                case 'LiF'
                    Current_FC = [pad(' ',52) '&' pad(' ',52)];
                case 'LiCl'
                    Current_Theory = [Current_Theory '~\cite{bach2009synthesis}'];
                    Current_FC = [pad(' ($\frac{1}{3}$, $\frac{2}{3}$, 0) ',52,'right') '&' pad(' ($\frac{1}{3}$, $\frac{2}{3}$, 0.379(1))',52,'right')];
                case 'LiBr'
                    Current_Theory = [Current_Theory '~\cite{liebold2008experimental}'];
                    Current_FC = [pad(' ($\frac{1}{3}$, $\frac{2}{3}$, 0) ',52,'right') '&' pad(' ($\frac{1}{3}$, $\frac{2}{3}$, 0.379(1))',52,'right')];
                case 'LiI'
                    Current_Theory = [Current_Theory '~\cite{fischer2004existiert}'];
                    Current_FC = [pad(' ',52) '&' pad(' ',52)];
                case 'NaCl'
                    Current_FC = [pad(' ',52) '&' pad(' ',52)];
            end
        elseif strcmp(Current_Structure,'Rocksalt') && strcmp(Current_Theory,'Experiment')
            Err = round(ExpData.(Current_Salt).da,1,'significant');
            n = 0;
            while (floor(Err*10^n)~=Err*10^n)
                n = n+1;
            end
            Errtxt = num2str(Err,'%1.0e');
            Current_a = [num2str(ExpData.(Current_Salt).a,['%5.' num2str(n) 'f']) '(' Errtxt(1) ')'];
            
            Current_Theory = [Current_Theory '~\cite{sirdeshmukh2013alkali}'];
            Current_FC = [pad(' (0, 0, 0) ',52,'right') '&' pad(' ($\frac{1}{2}$, $\frac{1}{2}$, $\frac{1}{2}$)',52,'right')];
        end
        
        % Check if structure has changed
        if ~strcmp(Current_Structure,Previous_Structure)
            
            % Previous Parameter Table Complete
            if Idx ~= 1
                                
                % Finish up all rows
                for i = 1:length(LinesText)
                    if i ~= length(LinesText)
                        LinesText{i} = [LinesText{i} '\\'];
                    end
                end

                % Build table
                z = z + 1;
                Parameter_Table{y,z} = [Line1 newline ...
                    Line2 newline ...
                    Line2b newline ...
                    LineAdj Adjusts newline ...
                    LineTitle];
                
                for i = 1:length(LinesText)
                    Parameter_Table{y,z} = [Parameter_Table{y,z} newline ...
                        LinesText{i}];
                end
                
                Parameter_Table{y,z} = [Parameter_Table{y,z} newline ...
                    Line5 newline ...
                    Line7 newline ...
                    Line8];
            end
            x = 0; % keep track of theory number
            
            Structure_Out = strrep(Current_Structure,'BetaBeO','$\beta$-BeO');
            Structure_Out = strrep(Structure_Out,'FiveFive','5-5');
            Structure_Out = strrep(Structure_Out,'Rocksalt','rocksalt');
            Structure_Out = strrep(Structure_Out,'Wurtzite','wurtzite');
            Structure_Out = strrep(Structure_Out,'Sphalerite','sphalerite');
            
            % Start new table
            Line1 = ['%Table: ' Current_Salt ' ' Current_Structure ' Structure Parameters ' Opt_Type];
            Line2 = '\begin{table}[H]';
            Line2b = '    \centering';
            LineAdj = '    \begin{tabular}';
            Line5 = '    \end{tabular}';
            Line7 = ['    \caption{\label{tab:' Current_Salt '_' Current_Structure '_' Label '} Optimized lattice parameters of ' Current_Salt ' in the ' Structure_Out ' crystal structure.}'];
            Line8 = '\end{table}';
            
            switch Current_Structure
                case 'Rocksalt'
                    Adjusts = '{l|lcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b = c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'Wurtzite'
                    Adjusts = '{l|llcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b}}',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'BetaBeO'
                    Adjusts = '{l|llcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b}}',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'CsCl'
                    Adjusts = '{l|lcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b = c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'FiveFive'
                    Adjusts = '{l|lllcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a}}',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{b}}',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'NiAs'
                    Adjusts = '{l|llcc}';
                    LineTitle = ['        ' pad(' ',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b}}',51,'right') '& ' pad('\multicolumn{1}{c}{\textbf{c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
                case 'Sphalerite'
                    Adjusts = '{l|lcc}';
                    LineTitle = ['        ' pad(' ',52,'right') '& ' pad('\multicolumn{1}{c}{\textbf{a = b = c}}',51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Metal ' Site}}'],51,'right') '& ' pad(['\multicolumn{1}{c}{\textbf{' Halide ' Site}}'],51,'right') '\\ \hline'];
            end
            LinesText = {};
            
        end
        x = x+1;
        
        if ~strcmp(Current_Salt,Previous_Salt)
            y = y + 1; % incriment salt tracker
            z = 0; % Reset structure tracker
        end
        
        if Remove_D3gCP
            Theory_Out = regexprep(Current_Theory,'-D3(-gCP){0,1}','');
        else
            Theory_Out = regexprep(Current_Theory,'(-D3)(-gCP){0,1}','-D3(BJ)$2');
        end
        Theory_Out = regexprep(Theory_Out,'JC(?!.)','JC (SPC/E)');
        Theory_Out = strrep(Theory_Out,'JC3P','JC (TIP3P)');
        Theory_Out = strrep(Theory_Out,'JC4P','JC (TIP4P$_{\text{EW}}$)');
        
        % Add in data
        LinesText{x} = ['        ' pad(['\textbf{' Theory_Out '}'],51,'right')];
        LinesText{x} = [LinesText{x} '& ' pad([Current_a Note],51,'right')];
        if strcmp(Current_Structure,'FiveFive')
            LinesText{x} = [LinesText{x} '& ' pad([Current_b Note],51,'right')];
        end
        if contained_in_cell(Current_Structure,{'Wurtzite' 'BetaBeO' 'FiveFive' 'NiAs'})
            LinesText{x} = [LinesText{x} '& ' pad([Current_c Note],51,'right')];
        end
        LinesText{x} = [LinesText{x} '&' Current_FC];
        
        Previous_Salt = Current_Salt;
        Previous_Structure = Current_Structure;
    end

    % Write last sturcture table
    % Finish up all rows
    for i = 1:length(LinesText)
        if i ~= length(LinesText)
            LinesText{i} = [LinesText{i} '\\'];
        end
    end

    % Build table
    z = z + 1;
    Parameter_Table{y,z} = [Line1 newline ...
        Line2 newline ...
        Line2b newline ...
        LineAdj Adjusts newline ...
        LineTitle];

    for i = 1:length(LinesText)
        Parameter_Table{y,z} = [Parameter_Table{y,end} newline ...
            LinesText{i}];
    end

    Parameter_Table{y,z} = [Parameter_Table{y,end} newline ...
        Line5 newline ...
        Line7 newline ...
        Line8];
    
    % Print Tables to file
    Text_Out = '';
    for i = 1:length(Energy_Table)
        Text_Out = [Text_Out Energy_Table{i} newline newline];
        
        for j = 1:size(Parameter_Table,2)
            Text_Out = [Text_Out Parameter_Table{i,j} newline newline ];
        end
    end
    waitbar(1,LoadbarObj,'Saving Tables to Disk...');
    % Save to file
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\n',Text_Out);
    fclose(fid);
    
    
    function Output_Structure = Shorten_Structure_Name(Input_Structure)
        Output_Structure = strrep(Input_Structure,'Rocksalt','Rock.');
        Output_Structure = strrep(Output_Structure,'Wurtzite','Wurtz.');
        Output_Structure = strrep(Output_Structure,'BetaBeO','$\beta$-BeO');
        Output_Structure = strrep(Output_Structure,'FiveFive','5-5');
        Output_Structure = strrep(Output_Structure,'Sphalerite','Sphal.');
    end
end