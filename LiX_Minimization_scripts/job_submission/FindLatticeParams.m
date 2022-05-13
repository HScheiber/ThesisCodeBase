function [a,b,c,FCmet,FChal] = FindLatticeParams(Salt,Structure,Model,home)

% How many decimal places of accuracy?
tol = 2;

% Load lattice energies quantum data
try
    X = load(fullfile(home,'data','CRYSTAL',...
    	[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
    QuantumData = X.Data;
    clearvars X
catch
    QuantumData = [];
end
    

% Load known empirical energies
try
    X = load(fullfile(home,'data','GROMACS',...
      [Salt '_' Structure '_Lattice_Energies.mat']),'Data');
    Data = X.Data;
    clearvars X
    
    EmpiricalData = Data.(Salt).(Structure).(Model);
    a_emp = [EmpiricalData{:,1}];
    b_emp = [EmpiricalData{:,2}];
    c_emp = [EmpiricalData{:,3}];
    FC_emp = {EmpiricalData{:,6}};
catch
    EmpiricalData = [];
    a_emp = [];
    b_emp = [];
    c_emp = [];
    FC_emp = {};
end

% Fields of quantum data available
if isempty(QuantumData)
    QuantumLabels = {};
else
    QuantumLabels = fields(QuantumData);
end


% Data of interest
Labend = Label_replace(Structure);
CLabel = [Salt Labend];

a_QM = [];
b_QM = [];
c_QM = [];
FC_QM = {};
% No matches found: return empty handed
if isempty(QuantumLabels)
    return
end

% Extract a, b, c, and FC sets from quantum mechanical data
for i=1:length(QuantumLabels)
    Label = QuantumLabels{i};

    DataOfInterest = QuantumData.(Label);
    
    % Match Basis set in data
    NBasis = size(DataOfInterest,1);
    for j=1:NBasis
        CheckingBasis = DataOfInterest{j,1};

        % Match basis set
        if strcmp(CheckingBasis,'pob-TZVP')
            a_QM = [a_QM DataOfInterest{j,2}{2}];
            b_QM = [b_QM DataOfInterest{j,2}{3}];
            c_QM = [c_QM DataOfInterest{j,2}{4}];
            FC_QM = [FC_QM DataOfInterest{j,2}{7}];
        else
            continue
        end
    end
end
N = length(a_QM);

% Find unique set of a, b, c, FC from quantum data
obj_array = cell(1,N);
W = ['%5.' num2str(tol) 'f'];
for i=1:N
    obj_array{i} = [num2str(a_QM(i),W) num2str(b_QM(i),W) num2str(c_QM(i),W) ...
        FC_QM{i}{1,1} num2str(mod(FC_QM{i}{1,2},1),W) num2str(mod(FC_QM{i}{1,3},1),W) num2str(mod(FC_QM{i}{1,4},1),W) ...
        FC_QM{i}{2,1} num2str(mod(FC_QM{i}{2,2},1),W) num2str(mod(FC_QM{i}{2,3},1),W) num2str(mod(FC_QM{i}{2,4},1),W)];
end
[~,ind,~] = unique(obj_array);
a_QM = a_QM(ind);
b_QM = b_QM(ind);
c_QM = c_QM(ind);
FC_QM = FC_QM(ind);

% Find un-matched data
Nemp = length(a_emp);
NQM = length(a_QM);

Emp_Obj = cell(1,Nemp);
QM_Obj = cell(1,NQM);
for i = 1:NQM
    QM_Obj{i} = [num2str(a_QM(i),W) num2str(b_QM(i),W) num2str(c_QM(i),W) ...
        FC_QM{i}{1,1} num2str(mod(FC_QM{i}{1,2},1),W) num2str(mod(FC_QM{i}{1,3},1),W) num2str(mod(FC_QM{i}{1,4},1),W) ...
        FC_QM{i}{2,1} num2str(mod(FC_QM{i}{2,2},1),W) num2str(mod(FC_QM{i}{2,3},1),W) num2str(mod(FC_QM{i}{2,4},1),W)];
end

for j = 1:Nemp
    Emp_Obj{j} = [num2str(a_emp(j),W) num2str(b_emp(j),W) num2str(c_emp(j),W) ...
        upper(FC_emp{j}{1,1}) num2str(mod(FC_emp{j}{1,2},1),W) num2str(mod(FC_emp{j}{1,3},1),W) num2str(mod(FC_emp{j}{1,4},1),W) ...
        upper(FC_emp{j}{2,1}) num2str(mod(FC_emp{j}{2,2},1),W) num2str(mod(FC_emp{j}{2,3},1),W) num2str(mod(FC_emp{j}{2,4},1),W)];
end
% QM objects that are not already in empirical data space
[~,indx] = setdiff(QM_Obj,Emp_Obj);

% Output lattice params
a = a_QM(indx);
b = b_QM(indx);
c = c_QM(indx);
FC = FC_QM(indx);

% Generate fractional coordinates in primary unit cell for output
M = length(FC);
FCmet = cell(M,1);
FChal = cell(M,1);
for i = 1:length(FC)
    [FCmet{i},FChal{i}] = UnitCell_FractionalCoords(...
        [FC{1}{1,2} FC{1}{1,3} FC{1}{1,4}],...
        [FC{1}{2,2} FC{1}{2,3} FC{1}{2,4}],Structure);
end

end


