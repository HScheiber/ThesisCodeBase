function Basis_Out = Ghost_Basis_BetaBeO(Salt,Met_Bas,Hal_Bas,BSSE_id)

[Metal,Halide] = Separate_Metal_Halide(Salt);
Metal_Z = elements('Symbol',Metal,'atomic_number');
Halide_Z = elements('Symbol',Halide,'atomic_number');

Met_ECP_txt = regexp(Met_Bas,'INPUT\n(([0-9]|\.)+ *)+','tokens','ONCE');
if isempty(Met_ECP_txt)
    Met_ghost = regexprep(Met_Bas,'^[0-9]+','100');
else
    Met_ECP_num = textscan(strtrim(Met_ECP_txt{1}),'%f');
    Met_ECP_Layers = sum(Met_ECP_num{1}(2:end));

    Met_ghost = regexprep(Met_Bas,['INPUT\n([0-9]| |\.)+?\n(.+?\n){' num2str(Met_ECP_Layers) '}'],'');
    Met_ghost = regexprep(Met_ghost,'^[0-9]+','100');
end

% Remove electrons
Met_ghost = regexprep(Met_ghost,'((0|\.)+) (([0-5]|\.)+) (([0-9]|\.)+) ([0-9]|\.)+ (([0-1]|\.)+)','$1 $2 $3 0.0 $4');

Hal_ECP_txt = regexp(Hal_Bas,'INPUT\n(([0-9]|\.)+ *)+','tokens','ONCE');
if isempty(Hal_ECP_txt)
    Hal_ghost = regexprep(Hal_Bas,'^[0-9]+','200');
else
    Hal_ECP_num = textscan(strtrim(Hal_ECP_txt{1}),'%f');
    Hal_ECP_Layers = sum(Hal_ECP_num{1}(2:end));

    Hal_ghost = regexprep(Hal_Bas,['INPUT\n([0-9]| |\.)+?\n(.+?\n){' num2str(Hal_ECP_Layers) '}'],'');
    Hal_ghost = regexprep(Hal_ghost,'^[0-9]+','200');
end

% Remove electrons
Hal_ghost = regexprep(Hal_ghost,'((0|\.)+) (([0-5]|\.)+) (([0-9]|\.)+) ([0-9]|\.)+ (([0-1]|\.)+)','$1 $2 $3 0.0 $4');

if strcmpi(BSSE_id,'Metal')
    Basis_Out = [Met_Bas newline Met_ghost newline Hal_ghost];
elseif strcmpi(BSSE_id,'Halide')
	Basis_Out = [Hal_Bas newline Met_ghost newline Hal_ghost];
end

end
