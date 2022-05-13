function [a,b,c,Coordinates] = FindCrystalParams(Salt,Crystal,Theory,Basis,...
    home,BSSE_id)

% How many decimal places of accuracy?
tol = 2;

% Get metal and halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

if strcmp(BSSE_id,'Metal')
    ion = Metal;
elseif strcmp(BSSE_id,'Halide')
    ion = Halide;
else
    error(['Unknown BSSE_id: ' BSSE_id])
end

% Load total energies data
X = load([home filesep 'data' filesep 'CRYSTAL' filesep ...
  Salt '_' Crystal '_Total_Energies.mat'],'Data');
Data = X.Data;
clearvars X

% Load known BSSE energies
X = load([home filesep 'data' filesep 'CRYSTAL' filesep ...
  ion '_BSSE_Energies.mat'],[BSSE_id '_BSSE_Data']);
BSSEData = X.([BSSE_id '_BSSE_Data']);
clearvars X

% Fields of data available
if isempty(Data)
    DataFields = {};
else
    DataFields = fields(Data);
end
if isempty(BSSEData)
    BSSE_DataFields = {};
else
    BSSE_DataFields = fields(BSSEData);
end

% Data of interest
Labend = Label_replace(Crystal);
CLabel = [Salt Labend '_' Theory];

% Search through fields for matches
matchedfields = regexp(DataFields,[CLabel '(_D[2-3].*)*$']);
matchedBSSEfields = regexp(BSSE_DataFields,CLabel);

matchedfieldslogical = ~cellfun(@isempty,matchedfields);
matchedBSSEfieldslogical = ~cellfun(@isempty,matchedBSSEfields);

Labels = DataFields(matchedfieldslogical);
BSSELabels = BSSE_DataFields(matchedBSSEfieldslogical);

a = [];
b = [];
c = [];
Coordinates = {};
% No match found: return empty handed
if isempty(Labels)
    return
end

for i=1:length(Labels)
    Label = Labels{i};
    Label_no_Disp = regexprep(Label,'(_D[2-3].*)$','');

    DataOfInterest = Data.(Label);
    if isfield(BSSEData,Label_no_Disp)
        BSSEDataOfInterest = BSSEData.(Label_no_Disp);
    else
        BSSEDataOfInterest = [];
    end
    
    % Match Basis set in data
    NBasis = size(DataOfInterest,1);
    matched = false;
    for j=1:NBasis
        CheckingBasis = DataOfInterest{j,1};

        % Match basis set
        if strcmp(Crystal,'Pair')
            if strcmp(CheckingBasis,Basis)
                matched = true;
                CellData = DataOfInterest{j,2};
                break
            end
        else
            if strcmp(CheckingBasis,'pob-TZVP') && strcmp(Basis,'def2-TZVPD')
                matched = true;
                CellData = DataOfInterest{j,2};
                break
            elseif strcmp(CheckingBasis,'pob-TZVP') && strcmp(Basis,'pob-TZVP')
                matched = false;
            elseif strcmp(CheckingBasis,'def2-TZVPD') && strcmp(Basis,'def2-TZVPD')
                matched = false;
            elseif strcmp(CheckingBasis,Basis)
                matched = true;
                CellData = DataOfInterest{j,2};
                break
            end
        end
    end
    
    % Match Basis set in already available BSSE data
    BSSEmatched = false;
    if ~isempty(BSSEDataOfInterest)
        NBasisBSSE = size(BSSEDataOfInterest,1);
        for j=1:NBasisBSSE
            CheckingBasis = BSSEDataOfInterest{j,1};
            
            if strcmp(CheckingBasis,Basis)
                BSSEmatched = true;
                BSSECellData = BSSEDataOfInterest{j,2};
                break
            end
        end
    end
    
    % Remove data that already exists
    if matched && BSSEmatched
        if strcmp(Crystal,'Pair')
            CellData{1,5} = setdiff(round(CellData{1,5}.*10^tol)./10^tol,round(BSSECellData{1,5}*10^tol)/10^tol);
        else
            [CellData{1,2},ia] = setdiff(round(CellData{1,2}.*10^tol)./10^tol,round(BSSECellData{1,2}*10^tol)/10^tol);
            CellData{1,3} = CellData{1,3}(ia);
            CellData{1,4} = CellData{1,4}(ia);
            CellData{1,7} = CellData{1,7}(ia);
        end
    end

    if matched
        if strcmp(Crystal,'Pair')
            a_current = CellData{1,5};
            b_current = [];
            c_current = [];
            Coord_current = {};
        else
            a_current = CellData{1,2};
            b_current = CellData{1,3};
            c_current = CellData{1,4};
            Coord_current = CellData{1,7};
        end
    else
        a_current = [];
        b_current = [];
        c_current = [];
        Coord_current = {};
    end
    
    a = [a a_current];
    b = [b b_current];
    c = [c c_current];
    Coordinates = [Coordinates Coord_current];
    
    % Remove duplicates
    [a,ind] = unique(round(a.*10^tol)./10^tol);
    if ~strcmp(Crystal,'Pair')
        b = b(ind);
        c = c(ind);
        Coordinates = Coordinates(ind);
    end
end
end