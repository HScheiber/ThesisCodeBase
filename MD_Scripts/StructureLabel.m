function Label = StructureLabel(Structure)
    switch lower(Structure)
        case 'rocksalt'
            Label = 'R';
        case 'wurtzite'
            Label = 'W';
        case 'sphalerite'
            Label = 'S';
        case 'cscl'
            Label = 'C';
        case 'nias'
            Label = 'N';
        case 'antinias'
            Label = 'A';
        case 'betabeo'
            Label = 'B';
        case 'fivefive'
            Label = 'F';
        case 'liquid'
            Label = 'L';
        case 'previous'
            Label = 'Pr';
    end
end