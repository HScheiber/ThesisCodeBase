function Bond_length = Bond_length_calculation(varargin)
if nargin > 4
    error('Too many input arguments.')
elseif nargin == 4
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
    structure = varargin{4};
elseif nargin == 2
    a = varargin{1};
    structure = varargin{2};
else
    error('Not enough input arguments.')
end

if strcmp(structure,'Rocksalt')
    Bond_length = a.*0.5;
elseif strcmp(structure,'Wurtzite')
    Bond_length = a.*((3./8)^(1/2));
elseif strcmp(structure,'Sphalerite')
    Bond_length = a.*sqrt(3)./4;
elseif strcmp(structure,'CsCl')
    Bond_length = a.*sqrt(3)./2;
elseif strcmp(structure,'NiAs')
    Bond_length = sqrt((a.^2)./3 + (c.^2)./16);
elseif strcmp(structure,'BetaBeO')
    Bond_length = a.*0.354;
elseif strcmp(structure,'FiveFive')
    Bond_length = a./2;
else
    error('Unknown Structure Type.')
end


end