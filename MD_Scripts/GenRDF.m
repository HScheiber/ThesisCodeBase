% Binsize is the granularity, in angstroms.
% Endpoint is the final point of the gR plot, first point is always 0, in angstroms
function [RDF,CRDF] = GenRDF(app,Input,BinSize,EndPoint)

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
R = 0:BinSize:EndPoint;

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
    RDF = [];
    CRDF = [];
    return
end

DT = [SubData{:,9}];
DatInd = ismember(DT,Data_Type);

if sum(DatInd) == 0
    warning(['Missing Requested Data for: ' Salt ' ' Structure ' ' Model]);
    RDF = [];
    CRDF = [];
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
[MX_BL,MM_BL,XX_BL,~] = Gen_Stars(CryStruc,Structure,Metal,Halide,EndPoint);
        
RDF.MX = zeros(1,length(R));
RDF.MM = zeros(1,length(R));
RDF.XX = zeros(1,length(R));
RDF.MA = zeros(1,length(R));
RDF.XA = zeros(1,length(R));
RDF.AA = zeros(1,length(R));

MXi = 1;
MXi_max = size(MX_BL,1);
MMi = 1;
MMi_max = size(MM_BL,1);
XXi = 1;
XXi_max = size(XX_BL,1);

for i = 1:length(R)-1
    while (MXi <= MXi_max) && ( R(i) > MX_BL(MXi,1) ) && ( MX_BL(MXi,1) <= R(i+1) )
        RDF.MX(i) = RDF.MX(i) + MX_BL(MXi,2);
        RDF.MA(i) = RDF.MA(i) + MX_BL(MXi,2);
        RDF.XA(i) = RDF.XA(i) + MX_BL(MXi,2);
        RDF.AA(i) = RDF.AA(i) + MX_BL(MXi,2);
        MXi = MXi + 1;
    end
    
    while (MMi <= MMi_max) && ( R(i) > MM_BL(MMi,1) ) && ( MM_BL(MMi,1) <= R(i+1) )
        RDF.MM(i) = RDF.MM(i) + MM_BL(MMi,2);
        RDF.MA(i) = RDF.MA(i) + MM_BL(MMi,2);
        RDF.AA(i) = RDF.AA(i) + MM_BL(MMi,2);
        MMi = MMi + 1;
    end
    
    while (XXi <= XXi_max) && ( R(i) > XX_BL(XXi,1) ) && ( XX_BL(XXi,1) <= R(i+1) )
        RDF.XX(i) = RDF.XX(i) + XX_BL(XXi,2);
        RDF.XA(i) = RDF.XA(i) + XX_BL(XXi,2);
        RDF.AA(i) = RDF.AA(i) + XX_BL(XXi,2);
        XXi = XXi + 1;
    end
end

% CRDFs
CRDF.MX = cumsum(RDF.MX);
CRDF.MM = cumsum(RDF.MM);
CRDF.XX = cumsum(RDF.XX);
CRDF.MA = cumsum(RDF.MA);
CRDF.XA = cumsum(RDF.XA);
CRDF.AA = cumsum(RDF.AA);

RDF.R = R;
CRDF.R = R;
end
