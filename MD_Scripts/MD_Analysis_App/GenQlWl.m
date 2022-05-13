
% l is a list of spherical harmonic primary numbers > 0 (even numbers only, e.g. l = 2:2:10!)
function [Ql,Wl] = GenQlWl(app,Input,NeighbourSelected,lims)

%% Where is data stored
datadir = [app.home filesep 'data' filesep 'GROMACS'];

%% Scaling parameters for TF and JC models
Salt = Input.Salt;
Structure = Input.Structure; %{'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'}; Initial structure
Model = Input.Model; % Input model(s) to use: JC, JC3P, JC4P, TF		 
Damping = Input.Damping; % Input damping functions to use: 0 to 6 are defined
TF_Paramset = Input.TF_Paramset; % choose from 0-3
S_D = Input.Scale_Dispersion; % Works for both JC and TF
S_R = Input.Scale_Repulsion; % Works for both JC and TF
S_MMD = Input.Scale_MM_Dispersion; % Works for both JC and TF
S_XXD = Input.Scale_XX_Dispersion; % Works for both JC and TF
S_MXD = Input.Scale_MX_Dispersion; % Works for both JC and TF
S_E = Input.Scale_Epsilon; % Scale all Epsilon (affects JC only)
S_S = Input.Scale_Sigma; % Scale all Sigma (affects JC only)
S_A = Input.Scale_Alpha; % Scale the repulsive exponential parameter alpha (affects TF only)
GAdjust_MX = Input.GAdjust_MX; % Gaussian Adjustment to the potential
GAdjust_MM = Input.GAdjust_MM; % Gaussian Adjustment to the potential
GAdjust_XX = Input.GAdjust_XX; % Gaussian Adjustment to the potential
Data_Type = Input.Data_Type; % 1 = cell optimization, 2 = full optimization

[Metal,Halide] = Separate_Metal_Halide(Salt);

% Build up the model name
G_Adj = '';
if sum(GAdjust_MX(:,1)) ~= 0
    G_Adj = [G_Adj '_GMX'];
end

if sum(GAdjust_MM(:,1)) ~= 0
    G_Adj = [G_Adj '_GMM'];
end

if sum(GAdjust_XX(:,1)) ~= 0
    G_Adj = [G_Adj '_GXX'];
end

if TF_Paramset ~= 0 && strcmp(Model,'TF')
    Model = [Model num2str(TF_Paramset)]; %#ok<*AGROW>
end

if Damping ~= 0
    Model = [Model 'd' num2str(Damping)];
end

if ~ismembertol(1.0,S_D,1e-5)
    mtxt = strrep(num2str(S_D,'%10.5f'),'-','N');
    Model = [Model '_D' mtxt];
end
if ~ismembertol(1.0,S_R,1e-5)
    mtxt = strrep(num2str(S_R,'%10.5f'),'-','N');
    Model = [Model '_R' mtxt];
end
if ~ismembertol(1.0,S_E,1e-5)
    mtxt = strrep(num2str(S_E,'%10.5f'),'-','N');
    Model = [Model '_E' mtxt];
end
if ~ismembertol(1.0,S_S,1e-5)
    mtxt = strrep(num2str(S_S,'%10.5f'),'-','N');
    Model = [Model '_S' mtxt];
end
if ~ismembertol(1.0,S_MMD,1e-5)
    mtxt = strrep(num2str(S_MMD,'%10.5f'),'-','N');
    Model = [Model '_MMD' mtxt];
end
if ~ismembertol(1.0,S_XXD,1e-5)
    mtxt = strrep(num2str(S_XXD,'%10.5f'),'-','N');
    Model = [Model '_XXD' mtxt];
end
if ~ismembertol(1.0,S_MXD,1e-5)
    mtxt = strrep(num2str(S_MXD,'%10.5f'),'-','N');
    Model = [Model '_MXD' mtxt];
end

if ~ismembertol(1.0,S_A,1e-5)
    mtxt = strrep(num2str(S_A,'%10.5f'),'-','N');
    Model = [Model '_A' mtxt];
end

Model = [Model G_Adj];
     
try
    % Load data of interest
    load(fullfile(datadir,[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
    SubData = Data.(Salt).(Structure).(Model);
catch
    warning(['Missing Data for: ' Salt ' ' Structure ' ' Model]);
    Ql = [];
    Wl = [];
    return
end

DT = [SubData{:,9}];
DatInd = ismember(DT,Data_Type);

if sum(DatInd) == 0
    warning(['Missing Requested Data for: ' Salt ' ' Structure ' ' Model]);
    Ql = [];
    Wl = [];
    return
end
            
CryStruc.a = [SubData{DatInd,1}];
CryStruc.b = [SubData{DatInd,2}];
CryStruc.c = [SubData{DatInd,3}];
CellAngles = SubData(DatInd,8);
CryStruc.FC = SubData(DatInd,6);   
if sum(DatInd) > 1
    En = [SubData{DatInd,7}];
    [~,midx] = min(En);
    CryStruc.a = CryStruc.a(midx);
    CryStruc.b = CryStruc.b(midx);
    CryStruc.c = CryStruc.c(midx);
    CryStruc.alpha = CellAngles{midx}(1);
    CryStruc.beta = CellAngles{midx}(2);
    CryStruc.gamma = CellAngles{midx}(3);
    CryStruc.FC = CryStruc.FC{midx};
else
    CryStruc.alpha = CellAngles{1}(1);
    CryStruc.beta = CellAngles{1}(2);
    CryStruc.gamma = CellAngles{1}(3);
    CryStruc.FC = CryStruc.FC{1};
end
       
% Generate coordination shells using structural parameters
if NeighbourSelected
    StartPoint = 0; % Angstroms
    EndPoint = 50; % Angstroms
    MaxNum = lims;

else
    StartPoint = lims(1); % Angstroms
    EndPoint = lims(2); % Angstroms
    MaxNum = Inf;
end
[MA_Clus,MX_Clus,MM_Clus,XA_Clus,XM_Clus,XX_Clus] = ...
    Gen_Stars_Cartesian(CryStruc,Structure,Metal,Halide,StartPoint,EndPoint,MaxNum);
        
l = app.l;
Ql.l = l;
Wl.l = l;

[Ql.MM,Wl.MM] = CalculateQlWl(MM_Clus,l);
[Ql.MX,Wl.MX] = CalculateQlWl(MX_Clus,l);
[Ql.MA,Wl.MA] = CalculateQlWl(MA_Clus,l);

[Ql.XX,Wl.XX] = CalculateQlWl(XX_Clus,l);
[Ql.XM,Wl.XM] = CalculateQlWl(XM_Clus,l);
[Ql.XA,Wl.XA] = CalculateQlWl(XA_Clus,l);

Ql.AA = (Ql.MA + Ql.XA)./2;
Wl.AA = (Wl.MA + Wl.XA)./2;

end
