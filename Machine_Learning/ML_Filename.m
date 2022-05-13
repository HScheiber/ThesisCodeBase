function Filename = ML_Filename(Salt,Structures,Model,Fix_Charge,Additivity,Loss_Options,Trial_ID)

Filename = [Salt '_' Model];

if Fix_Charge
    Filename = [Filename 'q'];
end

if strcmp(Model,'JC') && Additivity
    Filename = [Filename 'm'];
end

% Optional: add model ID
if ~isempty(Trial_ID)
    Filename = [Filename '_Model_' Trial_ID];
end

% Fill in the loss part
Filename = [Filename '_Loss'];
AE = '';
RE = '';
LP = '';
for idx = 1:length(Structures)
    Structure = Structures{idx};

    % Set the default weighting to 0
    if Loss_Options.(Structure).LE > 0
        AE = [AE Structure(1)];
    end
    if Loss_Options.(Structure).RLE > 0
        RE = [RE Structure(1)];
    end
    if (Loss_Options.(Structure).a + Loss_Options.(Structure).b + Loss_Options.(Structure).c) > 0
        LP = [LP Structure(1)];
    end
end

if ~isempty(AE)
    Filename = [Filename '_AE-' AE];
end
if ~isempty(RE)
    Filename = [Filename '_RE-' RE];
end
if ~isempty(LP)
    Filename = [Filename '_LP-' LP];
end

end