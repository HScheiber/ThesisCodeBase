function [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_STO_3G(home,current_crystal_type,Initialize_As_Ions)

load([home filesep 'basis_sets' filesep 'STO_3G.mat']) %#ok<LOAD>

% Modify Charges if needed
if strcmp(current_crystal_type,'Ion')
    Metal_Charge = Metal_Charge-1;
    Halide_Charge = Halide_Charge+1;
elseif strcmp(current_crystal_type,'Atom')
    % Do not modify charges
elseif Initialize_As_Ions
    Metal_Charge = Metal_Charge-1; %#ok<*NODEF>
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

end