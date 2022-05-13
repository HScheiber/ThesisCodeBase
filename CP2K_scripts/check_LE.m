function stop = check_LE(optimValues,~,Structure,M_En,X_En)

stop = false;

E_Conv = 2625.4996394799; % Hartree -> kJ/mol conversion factor
Total_energy = optimValues.fval;

switch Structure
    case {'Rocksalt' 'Sphalerite'}
        a = optimValues.x(1)*sqrt(2);
        b = optimValues.x(1)*sqrt(2);
        c = optimValues.x(1)*sqrt(2);
    case 'CsCl'
        a = optimValues.x(1);
        b = optimValues.x(1);
        c = optimValues.x(1);
    case 'FiveFive'
        a = optimValues.x(2);
        b = (3/sqrt(3))*optimValues.x(1);
        c = optimValues.x(1);
    case {'Wurtzite' 'NiAs' 'AntiNiAs'}
        a = optimValues.x(1);
        b = optimValues.x(1);
        c = optimValues.x(2);
    case 'BetaBeO'
        a = optimValues.x(1);
        b = optimValues.x(1);
        c = optimValues.x(2);
end

Lattice_Energy = (Total_energy - M_En - X_En)*E_Conv;

% Print
disp(repelem('*',60))
disp(['Iteration:' repelem(char(9),6) num2str(optimValues.iteration,'%i')])
disp(['Total Function Evaluations:' repelem(char(9),2) num2str(optimValues.funccount,'%i')])
if ~isempty(strtrim(optimValues.method))
    disp(['Update Method:' repelem(char(9),5) num2str(optimValues.method,'%i')])
end
disp(['Current Total Energy:' repelem(char(9),3) num2str(Total_energy,'%.8f') ' a.u.'])
disp(['Current Mesh Size:' repelem(char(9),4) num2str(optimValues.meshsize,'%.3e') ' ' char(0197)])
if ~isnan(optimValues.TolFun)
    disp(['Change in function value:' repelem(char(9),2) num2str(optimValues.TolFun,'%.3e') ' a.u.'])
end
if ~isnan(optimValues.TolX)
    disp(['Norm of change in x:' repelem(char(9),3) num2str(optimValues.TolX,'%.3e') ' ' char(0197)])
end
disp(['Current Lattice Energy:' repelem(char(9),3) 'E = ' num2str(Lattice_Energy,'%.6f') ' kJ/mol'])
disp(['Current Lattice Parameters:' repelem(char(9),2) 'a = ' num2str(a,'%.6f') ' ' char(0197) ...
    newline repelem(char(9),8) 'b = ' num2str(b,'%.6f') ' ' char(0197) ...
    newline repelem(char(9),8) 'c = ' num2str(c,'%.6f') ' ' char(0197)])
disp(repelem('*',60))

end