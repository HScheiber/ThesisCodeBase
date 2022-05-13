function Data = find_thermal_energies(jobtxt,div)

    % All energies are output in kJ/mol
    % Ouputs a structure with the following properties
    % Note that TS properties are actually -TS
    % Data.Pot_E
    % Data.ZPVE
    % Data.Thermal_E
    % Data.Thermal_G
    % Data.Total_G
    % Data.Trans_E
    % Data.Trans_TS
    % Data.Rot_E
    % Data.Rot_TS
    % Data.Vib_E
    % Data.Vib_TS
    % Data.Tot_E
    % Data.Tot_TS


    thermtxt = regexp(jobtxt,'Zero-point correction=(.*?\n){16}','match','once');
    
    if isempty(thermtxt)
        Data.Pot_E = nan;
        Data.ZPVE = nan;
        Data.Thermal_E = nan;
        Data.Thermal_G = nan;
        Data.Total_G = nan;
        Data.Trans_E = nan;
        Data.Trans_TS = nan;
        Data.Rot_E = nan;
        Data.Rot_TS = nan;
        Data.Vib_E = nan;
        Data.Vib_TS = nan;
        Data.Tot_E = nan;
        Data.Tot_TS = nan;
        return
    end
    
    if nargin < 2
        div = 1;
    end
    kJ_mol_per_Ha = 2625.4996394799;
    J_cal = 4.184;
    T = 298.15;

    ZPVE = regexp(thermtxt,'Zero-point correction= +([0-9]|\.|-)+','tokens','once');
    Data.ZPVE = str2double(ZPVE{1})*kJ_mol_per_Ha/div; % kJ/mol
    
    Thermal_E = regexp(thermtxt,'Thermal correction to Energy= +([0-9]|\.|-)+','tokens','once');
    Data.Thermal_E = str2double(Thermal_E{1})*kJ_mol_per_Ha/div;
    
    Thermal_G = regexp(thermtxt,'Thermal correction to Gibbs Free Energy= +([0-9]|\.|-)+','tokens','once');
    Data.Thermal_G = str2double(Thermal_G{1})*kJ_mol_per_Ha/div;
    
    Total_G = regexp(thermtxt,'Sum of electronic and thermal Free Energies= +([0-9]|\.|-)+','tokens','once');
    Data.Total_G = str2double(Total_G{1})*kJ_mol_per_Ha/div;
    
    Sum_ZPVE = regexp(thermtxt,'Sum of electronic and zero-point Energies= +([0-9]|\.|-)+','tokens','once');
    Sum_ZPVE = str2double(Sum_ZPVE{1})*kJ_mol_per_Ha/div;
    Data.Pot_E = Sum_ZPVE - Data.ZPVE;
    
    Trans_ES = regexp(thermtxt,'Translational +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
    Data.Trans_E = str2double(Trans_ES{1})*J_cal/div; % units in kJ/mol
    Data.Trans_TS = (-str2double(Trans_ES{3})*J_cal*T/1000)/div; % units in Jk/mol-K
    
    Rot_ES = regexp(thermtxt,'Rotational +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
    Data.Rot_E = str2double(Rot_ES{1})*J_cal/div; % units in kJ/mol
    Data.Rot_TS = (-str2double(Rot_ES{3})*J_cal*T/1000)/div; % units in Jk/mol-K

    Vib_ES = regexp(thermtxt,'Vibrational +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
    Data.Vib_E = str2double(Vib_ES{1})*J_cal/div; % units in kJ/mol
    Data.Vib_TS = (-str2double(Vib_ES{3})*J_cal*T/1000)/div; % units in kJ/mol-K

    Tot_ES = regexp(thermtxt,'Total +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
    Data.Tot_E = str2double(Tot_ES{1})*J_cal/div; % units in kJ/mol
    Data.Tot_TS = (-str2double(Tot_ES{3})*J_cal*T/1000)/div; % units in kJ/mol

end