function [Li,Na,F,Cl,Br] = Prune_Scale_cc_pV5Z(home,current_crystal_type,BSSE_id,Initialize_As_Ions)

load([home filesep 'basis_sets' filesep 'cc_pV5Z_CRYSTAL17.mat']) %#ok<LOAD>
RegS = ' ([1-9]|10)\.0 '; % Regex string to replace charges

% Modify Charges if needed
if strcmp(current_crystal_type,'Ion')
    Metal_Charge = Metal_Charge-1; %#ok<*NODEF>
    Halide_Charge = Halide_Charge+1;
elseif strcmp(current_crystal_type,'Atom')
    % Do not modify charges
elseif Initialize_As_Ions
    Metal_Charge = Metal_Charge-1;
    Halide_Charge = Halide_Charge+1;
end

% Replace Charges back into text
Li = strrep(Li,'##CHG##',num2str(Metal_Charge,'%4.1f'));
Na = strrep(Na,'##CHG##',num2str(Metal_Charge,'%4.1f'));
F = strrep(F,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Cl = strrep(Cl,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Br = strrep(Br,'##CHG##',num2str(Halide_Charge,'%4.1f'));

% For custom pair BSSE
if strcmp(BSSE_id,'Halide') && strcmp(current_crystal_type,'Pair')
    Li_ATNUM = 0;
    Na_ATNUM = 0;
    
    % Remove electrons
    Li = regexprep(Li,RegS,' 0.0 ');
    Na = regexprep(Na,RegS,' 0.0 ');
elseif strcmp(BSSE_id,'Metal') && strcmp(current_crystal_type,'Pair')
    F_ATNUM = 0;
    Cl_ATNUM = 0;
    Br_ATNUM = 0;
    
    % Remove electrons
    F = regexprep(F,RegS,' 0.0 ');
    Cl = regexprep(Cl,RegS,' 0.0 ');
    Br = regexprep(Br,RegS,' 0.0 ');
end

Li = strrep(Li,'##ATNUM##',num2str(Li_ATNUM,'%u'));
Na = strrep(Na,'##ATNUM##',num2str(Na_ATNUM,'%u'));
F = strrep(F,'##ATNUM##',num2str(F_ATNUM,'%u'));
Cl = strrep(Cl,'##ATNUM##',num2str(Cl_ATNUM,'%u'));
Br = strrep(Br,'##ATNUM##',num2str(Br_ATNUM,'%u'));

end