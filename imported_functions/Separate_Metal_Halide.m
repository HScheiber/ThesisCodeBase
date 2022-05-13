function [Metal,Halide] = Separate_Metal_Halide(Salt)
% Preprocessing
[match,~] = regexp(Salt,{'M' 'Li' 'Na' 'K' 'Rb' 'Cs' 'Fe' 'Al' 'F' 'Cl' 'Br' 'I' 'At' 'X' 'O'},'match','tokens');
Matched_Metals = match(1:8);
Matched_Halides = match(9:end);
matches_metals = Matched_Metals(~cellfun('isempty',Matched_Metals));
matches_halides = Matched_Halides(~cellfun('isempty',Matched_Halides));
if ~isempty(matches_metals)
    Metal = matches_metals{end}{1};
else
    Metal = '';
end

if ~isempty(matches_halides)
    Halide = matches_halides{end}{1};
else
    Halide = '';
end
end