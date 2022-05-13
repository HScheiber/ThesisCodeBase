function Geometry = SymmetryAdapt(Geometry,Structure)
switch Structure
    case 'BetaBeO'
        %c = (1/sqrt(3)).*a_in;
        Geometry.b = Geometry.a;
        
    case 'CsCl'
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;

    case 'FiveFive'
        %b = (3/2).*a_in; % Theoretical b/a = 3/2
        %c = (sqrt(3)/2).*a_in; % Theoretical c/a = cosd(30)
        if ~Geometry.Conv
            Geometry.b = Geometry.a;
        end
        
    case {'NiAs' 'AntiNiAs'}
        %x = 1.391666666667;
        % x.*a_in; 1.7316 LiF LSDA, 1.7275 LiF PBE
        Geometry.b = Geometry.a;
    case 'Rocksalt'
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;
        
    case 'Sphalerite'
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;
        
    case 'Wurtzite'
        %c = sqrt(8/3).*a_in; % Perfect Wurtzite c/a = sqrt(8/3);
        Geometry.b = Geometry.a;
end

end