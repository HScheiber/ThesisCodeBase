% This function will center fractional coordinates of a crystal so that none are along the unit cell border
function Geometry = CenterCoordinates(Geometry,Structure)

    % Conventional structure
    if Geometry.Conv
        switch lower(Structure)
            case 'rocksalt'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.25 0.25 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.25 0.25 0.25],1);
            case 'wurtzite'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.3125],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.3125],1);
            case 'sphalerite'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.125 0.125 0.125],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.125 0.125 0.125],1);
            case 'betabeo'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.25],1);
            case {'nias' 'antinias'}
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [1/6 1/6 1/8],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [1/6 1/6 1/8],1);
            case 'fivefive'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.25],1);
            case 'cscl'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.25 0.25 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.25 0.25 0.25],1);
        end

    % Primitive structure
    else
        switch lower(Structure)
            case 'rocksalt'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.25 0.25 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.25 0.25 0.25],1);
            case 'wurtzite'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.3125],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.3125],1);
            case 'sphalerite'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.3750 0.3750 0.3750],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.3750 0.3750 0.3750],1);
            case 'betabeo'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.25],1);
            case {'nias' 'antinias'}
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [1/6 1/6 1/8],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [1/6 1/6 1/8],1);
            case 'fivefive'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0 0 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0 0 0.25],1);
            case 'cscl'
                Geometry.FC_Metal = mod(Geometry.FC_Metal + [0.25 0.25 0.25],1);
                Geometry.FC_Halide = mod(Geometry.FC_Halide + [0.25 0.25 0.25],1);
        end
    end
    Geometry.FC = [Geometry.FC_Metal; Geometry.FC_Halide];
end