
function View_Gmx_Table(Table_Loc,Q,C)


    nm_per_m = 1e+9; % nm per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1


    %% Read columns of data as text:
    % For more information, see the TEXTSCAN documentation.
    %formatSpec = '%16s%19s%17s%20s%18s%19s%s%[^\n\r]';
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

    %% Open the text file.
    fileID = fopen(Table_Loc,'r');

    %% Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

    %% Close the text file.
    fclose(fileID);

%     %% Convert the contents of columns containing numeric text to numbers.
%     % Replace non-numeric text with NaN.
%     raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
%     for col=1:length(dataArray)-1
%         raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
%     end
%     numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% 
%     for col=[1,2,3,4,5,6,7]
%         % Converts text in the input cell array to numbers. Replaced non-numeric text with NaN.
%         rawData = dataArray{col};
%         for row=1:size(rawData, 1)
%             % Create a regular expression to detect and remove non-numeric prefixes and suffixes.
%             regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%             try
%                 result = regexp(rawData(row), regexstr, 'names');
%                 numbers = result.numbers;
% 
%                 % Detected commas in non-thousand locations.
%                 invalidThousandsSeparator = false;
%                 if numbers.contains(',')
%                     thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
%                     if isempty(regexp(numbers, thousandsRegExp, 'once'))
%                         numbers = NaN;
%                         invalidThousandsSeparator = true;
%                     end
%                 end
%                 % Convert numeric text to numbers.
%                 if ~invalidThousandsSeparator
%                     numbers = textscan(char(strrep(numbers, ',', '')), '%f');
%                     numericData(row, col) = numbers{1};
%                     raw{row, col} = numbers{1};
%                 end
%             catch
%                 raw{row, col} = rawData{row};
%             end
%         end
%     end


    %% Replace non-numeric cells with NaN
    raw = dataArray;
%     R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
%     raw(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    U = table;
    U.r = cell2mat(raw(:, 1));
    U.f = cell2mat(raw(:, 2));
    U.df = cell2mat(raw(:, 3));
    U.h = cell2mat(raw(:, 4));
    U.dh = cell2mat(raw(:, 5));
    U.g = cell2mat(raw(:, 6));
    U.dg = cell2mat(raw(:, 7));

    %% Clear temporary variables
    clearvars filename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
    U_Total = k_0*(e_c^2).*Q.*U.f + C.*U.h + U.g;
    dU_Total = -k_0*(e_c^2).*Q.*U.df - C.*U.dh  - U.dg;
    %ref_U = %k_0*(e_c^2).*Q.*(1./U.r);
    
    hold on
    plot(U.r,U_Total,'k')
    hold on
    plot(U.r,dU_Total,'r')
    %hold on
    %plot(U.r,ref_U)
    ylim([-800 100])
    xlim([0 1])
end