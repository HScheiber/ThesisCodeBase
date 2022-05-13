function [Coords,Angles,ABC,Symmetry] = cp2k_Structures_fixed(Structure,Metal,Halide,Params)

% Possible Structures
% Rocksalt
% Wurtzite
% CsCl
% FiveFive
% BetaBeO
% Sphalerite
% NiAs
% AntiNiAs

switch lower(Structure)
    case 'rocksalt'
        Coords =   ['		' Metal '         0         0         0' newline ...
                    '		' Halide '   0.5000    0.5000    0.5000'];
        Angles = [60 60 60];
        Symmetry = 'RHOMBOHEDRAL';
        abc = [Params(1) Params(1) Params(1)];
    case 'wurtzite'
        Coords =   ['		' Metal '     0.666666666667 0.666666666667 0.375000000000' newline ...
                    '		' Metal '     0.333333333333 0.333333333333 0.875000000000' newline ...
                    '		' Halide '    0.666666666667 0.666666666667 0.000000000000' newline ...
                    '		' Halide '    0.333333333333 0.333333333333 0.500000000000'];
        Angles = [90 90 60];
        Symmetry = 'HEXAGONAL';
        abc = [Params(1) Params(1) Params(2)];
    case 'fivefive'
        Coords =   ['		' Metal '     0.666666666667 0.666666666667 0.500000000000' newline ...
                    '		' Metal '     0.333333333333 0.333333333333 0.000000000000' newline ...
                    '		' Halide '    0.666666666667 0.666666666667 0.000000000000' newline ...
                    '		' Halide '    0.333333333333 0.333333333333 0.500000000000'];
        Angles = [90 90 60];
        Symmetry = 'HEXAGONAL';
        abc = [Params(1) Params(1) Params(2)];
    case 'cscl'
        Coords =   ['		' Metal '         0         0         0' newline ...
                    '		' Halide '   0.5000    0.5000    0.5000'];
        Angles = [90 90 90];
        Symmetry = 'CUBIC';
        abc = [Params(1) Params(1) Params(1)];
    case 'betabeo'
        Coords =   ['		' Metal '    0.3360    0.6640    0.0000' newline ...
                    '		' Metal '    0.6640    0.3360    0.0000' newline ...
                    '		' Metal '    0.8360    0.8360    0.5000' newline ...
                    '		' Metal '    0.1640    0.1640    0.5000' newline ...
                    '		' Halide '   0.3100    0.3100    0.0000' newline ...
                    '		' Halide '   0.6900    0.6900    0.0000' newline ...
                    '		' Halide '   0.8100    0.1900    0.5000' newline ...
                    '		' Halide '   0.1900    0.8100    0.5000'];
        Angles = [90 90 90];
        Symmetry = 'TETRAGONAL';
        abc = [Params(1) Params(1) Params(2)];
    case 'sphalerite'
        Coords =   ['		' Metal '    0.0000    0.0000    0.0000' newline ...
                    '		' Halide '   0.2500    0.2500    0.2500'];
        Angles = [60 60 60];
        Symmetry = 'RHOMBOHEDRAL';
        abc = [Params(1) Params(1) Params(2)];
    case 'nias'
        Coords =   ['		' Metal '    0.0000            0.0000            0.0000' newline ...
                    '		' Metal '    0.0000            0.0000            0.5000' newline ...
                    '		' Halide '   0.333333333333    0.333333333333    0.2500' newline ...
                    '		' Halide '   0.666666666667    0.666666666667    0.7500'];
        Angles = [90 90 60];
        Symmetry = 'HEXAGONAL';
        abc = [Params(1) Params(1) Params(2)];
    case 'antinias'
        Coords =   ['		' Metal '   0.333333333333    0.333333333333    0.2500' newline ...
                    '		' Metal '   0.666666666667    0.666666666667    0.7500' newline ...
                    '		' Halide '    0.0000            0.0000            0.0000' newline ...
                    '		' Halide '    0.0000            0.0000            0.5000'];
        Angles = [90 90 60];
        Symmetry = 'HEXAGONAL';
        abc = [Params(1) Params(1) Params(2)];
end

ABC = num2str(abc,'%1.12f   ');
Angles = num2str(Angles,'%3i ');

end