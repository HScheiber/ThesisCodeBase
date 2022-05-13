function [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_def2_TZVPD(home,current_crystal_type,BSSE_id,Initialize_As_Ions)

load([home filesep 'basis_sets' filesep 'def2_TZVPD_CRYSTAL17.mat']) %#ok<LOAD>
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
K = strrep(K,'##CHG##',num2str(Metal_Charge,'%4.1f'));
Rb = strrep(Rb,'##CHG##',num2str(Metal_Charge,'%4.1f'));
F = strrep(F,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Cl = strrep(Cl,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Br = strrep(Br,'##CHG##',num2str(Halide_Charge,'%4.1f'));
I = strrep(I,'##CHG##',num2str(Halide_Charge,'%4.1f'));

% For custom pair BSSE
if strcmp(BSSE_id,'Halide') && strcmp(current_crystal_type,'Pair')
    Li_ATNUM = 0;
    Na_ATNUM = 0;
    K_ATNUM = 0;
    Rb_ATNUM = 0;
    Rb = strrep(Rb,['##INPUT##' newline],'');
    
    % Remove electrons
    Li = regexprep(Li,RegS,' 0.0 ');
    Na = regexprep(Na,RegS,' 0.0 ');
    K = regexprep(K,RegS,' 0.0 ');
    Rb = regexprep(Rb,RegS,' 0.0 ');
elseif strcmp(BSSE_id,'Metal') && strcmp(current_crystal_type,'Pair')
    F_ATNUM = 0;
    Cl_ATNUM = 0;
    Br_ATNUM = 0;
    I_ATNUM = 0;
    I = strrep(I,['##INPUT##' newline],'');
    
    % Remove electrons
    F = regexprep(F,RegS,' 0.0 ');
    Cl = regexprep(Cl,RegS,' 0.0 ');
    Br = regexprep(Br,RegS,' 0.0 ');
    I = regexprep(I,RegS,' 0.0 ');
end

Li = strrep(Li,'##ATNUM##',num2str(Li_ATNUM,'%u'));
Na = strrep(Na,'##ATNUM##',num2str(Na_ATNUM,'%u'));
K = strrep(K,'##ATNUM##',num2str(K_ATNUM,'%u'));
Rb = strrep(Rb,'##ATNUM##',num2str(Rb_ATNUM,'%u'));
Rb = strrep(Rb,'##INPUT##',Rb_INPUT);
F = strrep(F,'##ATNUM##',num2str(F_ATNUM,'%u'));
Cl = strrep(Cl,'##ATNUM##',num2str(Cl_ATNUM,'%u'));
Br = strrep(Br,'##ATNUM##',num2str(Br_ATNUM,'%u'));
I = strrep(I,'##ATNUM##',num2str(I_ATNUM,'%u'));
I = strrep(I,'##INPUT##',I_INPUT);

end