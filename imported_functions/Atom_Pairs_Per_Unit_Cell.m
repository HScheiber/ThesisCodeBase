function N = Atom_Pairs_Per_Unit_Cell(Structure)
if strcmp(Structure,'BetaBeO')
    N = 4;
elseif strcmp(Structure,'CsCl')
    N = 1;
elseif strcmp(Structure,'FiveFive')
    N = 4;
elseif strcmp(Structure,'NiAs')
    N = 2;
elseif strcmp(Structure,'Rocksalt')
    N = 4;
elseif strcmp(Structure,'Sphalerite')
    N = 4;
elseif strcmp(Structure,'Wurtzite')
    N = 2;
elseif strcmp(Structure,'Pair')
    N = 1;
else
    error(['Unknown Input Structure: ' Structure])
end