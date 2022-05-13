function N = Atoms_Per_Unit_Cell(Structure)
if strcmp(Structure,'BetaBeO')
    N = 8;
elseif strcmp(Structure,'CsCl')
    N = 2;
elseif strcmp(Structure,'FiveFive')
    N = 8;
elseif strcmp(Structure,'NiAs')
    N = 4;
elseif strcmp(Structure,'Rocksalt')
    N = 8;
elseif strcmp(Structure,'Sphalerite')
    N = 8;
elseif strcmp(Structure,'Wurtzite')
    N = 4;
elseif strcmp(Structure,'Pair')
    N = 2;
else
    error(['Unknown Input Structure: ' Structure])
end