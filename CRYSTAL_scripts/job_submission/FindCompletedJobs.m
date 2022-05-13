function SkipSwitch = FindCompletedJobs(Salt,Structure,Theory,Basis,...
    Dispersion,gCP,P1_Symmetry,External_P,Cell_optimization,...
    opt_type,opt_conv,QHA,Calc_Vib,Supercell,a,home)

if Cell_optimization
    switch opt_type
        case 'CELLONLY'
            DataType = 1;
        case {'FULLOPTG' 'ITATOCEL'}
            if P1_Symmetry
                DataType = 4;
            elseif External_P == 1
                DataType = 5;
            else
                DataType = 2;
            end
        case 'ATOMONLY'
            DataType = 3;
    end
else
    DataType = 0;
end

% Initialize the Skip switch as false (meaning do not skip this data)
SkipSwitch = false;

% Load data from disk
X = load([home filesep 'data' filesep 'CRYSTAL' filesep ...
  Salt '_' Structure '_Total_Energies.mat'],'Data');
Data = X.Data;
clearvars X

if ~isempty(Dispersion)
    Dispersion_lab = ['_' Dispersion];
else
    Dispersion_lab = '';
end

if ~isempty(gCP)
    gCP_lab = ['_' gCP];
else
    gCP_lab = '';
end

% Prepare label
Label = [Salt Label_replace(Structure) '_' Theory Dispersion_lab gCP_lab];

% Data of interest if available
if isfield(Data,Label)
    DOI = Data.(Label);
else
    % If data of interest does not exist at this level, we will not skip.
    return
end

% Next find basis set match
DOI_dat = [];
for i = 1:size(DOI,1)
    if strcmpi(DOI{i,1},Basis)
        DOI_dat = DOI{i,2};
        DOI_dat_Vib = DOI{i,5};
        DOI_dat_QHA = DOI{i,6};
        break
    end
end
if isempty(DOI_dat)
    % If data of interest does not exist at this level, we will not skip.
    return
end

% If run is vibrational anaylsis or QHA
if QHA
    if isempty(DOI_dat_QHA)
        % If data of interest does not exist at this level, we will not skip.
        return
    end
    % If data exists, look for matching supercell
    for i = 1:size(DOI_dat_QHA,1)
        if Supercell == DOI_dat_QHA{i,1}
            % Match found
            SkipSwitch = true;
            return
        end
    end
    % Otherwise, data is not found and I can return
    return
elseif Calc_Vib
    if isempty(DOI_dat_Vib)
        % If data of interest does not exist at this level, we will not skip.
        return
    end
    
    % If data exists, look for matching supercell
    for i = 1:size(DOI_dat_Vib,1)
        if Supercell == DOI_dat_Vib{i,1}
            % Match found
            SkipSwitch = true;
            return
        end
    end
    % Otherwise, data is not found and I can return
    return
end

% Back to geometry optimization or fixed cell
% Find only chosen data type
Conv_Data = cell(1,13);
Conv_Ind = (DOI_dat{1,10} == DataType);

% If no data of interest at this level, we will not skip
if sum(Conv_Ind) == 0
    return
end

% Gather the data of interest
for i = 2:14
    Conv_Data{i-1} = DOI_dat{1,i}(Conv_Ind);
end

if isempty(Conv_Data{1,1})
    % If data of interest does not exist at this level, we will not skip.
    return
end

% Match the a parameter
if DataType == 0
	SkipSwitch = ismembertol(a,[Conv_Data{:,1}],1e-5);
% Match the convergence criteria
elseif Cell_optimization
    Delta_toldee = abs(Conv_Data{13}(:) - opt_conv);
    
    if any(Delta_toldee < 1e-5) % Match within a tolerance
        SkipSwitch = true;
        return
    else
        return
    end
else % otherwise a match has been found, skip this one
    SkipSwitch = true;
end

end
