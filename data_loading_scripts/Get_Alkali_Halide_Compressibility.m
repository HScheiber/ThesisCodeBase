% Solid Alkali Halide Compressibility data from table 4: https://www-nature-com.ezproxy.library.ubc.ca/articles/srep03068/tables/4
% Isothermal Molten Salt compressibilities from: https://aip.scitation.org/doi/full/10.1063/1.4822097
function Compressibility = Get_Alkali_Halide_Compressibility(Salt,varargin) % Output should be in Bar^-1

p = inputParser;
p.FunctionName = 'Get_Alkali_Halide_Compressibility';
addOptional(p,'Isotropy','isotropic')
addOptional(p,'Molten',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
addOptional(p,'ScaleFactor',1)
parse(p,varargin{:});
Molten = p.Results.Molten;
Isotropy = p.Results.Isotropy;
SF = p.Results.ScaleFactor;

% Unit conversion
GPa_per_Bar = 1/10000;

% Liquid salts, isothermal compressibilities in units of GPa^(-1)
switch Salt
    case 'LiF'
        Liq_Compressibility = 0.093; 
    case 'LiCl'
        Liq_Compressibility = 0.216;
    case 'LiBr'
        Liq_Compressibility = 0.235;
    case 'LiI'
        Liq_Compressibility = 0.312;

    case 'NaF'
        Liq_Compressibility = 0.133;
    case 'NaCl'
        Liq_Compressibility = 0.343;
    case 'NaBr'
        Liq_Compressibility = 0.361;
    case 'NaI'
        Liq_Compressibility = 0.436;

    case 'KF'
        Liq_Compressibility = 0.186;
    case 'KCl'
        Liq_Compressibility = 0.442;
    case 'KBr'
        Liq_Compressibility = 0.465;
    case 'KI'
        Liq_Compressibility = 0.572;


    case 'RbF'
        Liq_Compressibility = 0.176;
    case 'RbCl'
        Liq_Compressibility = 0.429;
    case 'RbBr'
        Liq_Compressibility = 0.499;
    case 'RbI'
        Liq_Compressibility = 0.607;

    case 'CsF'
        Liq_Compressibility = 0.228;
    case 'CsCl'
        Liq_Compressibility = 0.461;
    case 'CsBr'
        Liq_Compressibility = 0.584;
    case 'CsI'
        Liq_Compressibility = 0.690;
    otherwise
        Liq_Compressibility = 0.4500;
        % H2O = 4.76e-5 bar^(-1), gromacs default = 4.5e-5 Bar^(-1)
end
Liq_Compressibility = Liq_Compressibility*GPa_per_Bar; % Bar^(-1)

% Bulk modulus of solid alkali halides in units of GPa
switch Salt
    case 'LiF'
        B = 76.4; 
    case 'LiCl'
        B = 32.9;
    case 'LiBr'
        B = 26.0;
    case 'LiI'
        B = 19.3;

    case 'NaF'
        B = 47.1;
    case 'NaCl'
        B = 23.9;
    case 'NaBr'
        B = 19.6;
    case 'NaI'
        B = 14.9;

    case 'KF'
        B = 28.7;
    case 'KCl'
        B = 16.5;
    case 'KBr'
        B = 13.9;
    case 'KI'
        B = 11.1;


    case 'RbF'
        B = 24.1;
    case 'RbCl'
        B = 14.0;
    case 'RbBr'
        B = 12.0;
    case 'RbI'
        B = 9.6;

    case 'CsF'
        B = 19.3;
    case 'CsCl'
        B = 17.6;
    case 'CsBr'
        B = 15.4;
    case 'CsI'
        B = 12.5;
    otherwise
        B = 1/0.4500; % Gromacs default bulk modulus
end
Sol_Compressibility = (1/B)*GPa_per_Bar; % Bar^(-1)

switch lower(Isotropy)
    % Single phase case
    case 'isotropic'
        if Molten
            Compressibility = Liq_Compressibility;
        else
            Compressibility = Sol_Compressibility;
        end
    % Solid-liquid case: Solid in first dimension and liquid in second dimension
    case 'semiisotropic'
        if Molten
            Compressibility = [1 1].*Liq_Compressibility;
        else
            Compressibility = [1 1].*Sol_Compressibility;
        end
    case 'anisotropic'
        if Molten
            Compressibility = [1 1 1 0 0 0].*Liq_Compressibility;
        else
            Compressibility = [1 1 1 0 0 0].*Sol_Compressibility;
        end
    case 'surface-tension'
        if Molten
            Compressibility = [1 1].*Liq_Compressibility;
        else
            Compressibility = [1 1].*Sol_Compressibility;
        end
end

Compressibility = SF.*Compressibility;
end