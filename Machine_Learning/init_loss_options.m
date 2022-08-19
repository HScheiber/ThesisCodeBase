function Loss_Options = init_loss_options

    % Possible options:
    % Absolue lattice energy of all structures
    % Relative lattice energy of all structures besides rocksalt
    % Lattice parameters of all structures
    % Weighting for each thing
    % L1 vs L2 regularization
    Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'FiveFive' 'Sphalerite' 'AntiNiAs' 'BetaBeO' 'CsCl'}; % 'BetaBeO' 'AntiNiAs' 'CsCl'

    Loss_Options.regularization = 'L2';
    
    % When true: use the experimental rocksalt lattice energy as RS LE
    % target, shift relative energies of all other structures to match the DFT gap
    % When false: use the DFT-calculated rocksalt lattice energy as RS LE
    % target, use DFT energies of all other structures.
    Loss_Options.Experimental_LE = true;
    
    % When true: use experimental lattice parameters for rocksalt and
    % wurtzite (where available)
    % When false: use DFT-calculated lattice parameters
    Loss_Options.Experimental_LP = true;
    
    % High temperature experimental properties of the liquid and solid
    Loss_Options.Fusion_Enthalpy  = 0; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
    Loss_Options.MP_Volume_Change = 0; % Fitting the experimental change in volume due to melting at the experimental MP
    Loss_Options.Liquid_MP_Volume = 0; % Fitting the experimental volume per formula unit at the experimental MP
    Loss_Options.Solid_MP_Volume  = 0; % Fitting the experimental volume of the experimental solid structure at the experimental MP
    Loss_Options.Liquid_DM_MP = 0; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
    Loss_Options.MP  = 0; % Fitting the experimental MP, using the experimental structure as the solid
    
    % Properties of the solid structures
    for idx = 1:length(Structures)
        Structure = Structures{idx};
        
        % Set the default weighting to 0
        Loss_Options.(Structure).LE  = 0;
        Loss_Options.(Structure).LE_freedom = 0; % allowable "give" in the lattice energy target, in kJ/mol
        Loss_Options.(Structure).RLE = 0; % weight
        Loss_Options.(Structure).RLE_freedom = 0; % allowable "give" in the relative lattice energy target, in kJ/mol
        Loss_Options.(Structure).a   = 0;
        Loss_Options.(Structure).b   = 0;
        Loss_Options.(Structure).c   = 0;
        Loss_Options.(Structure).V   = 0; % Volume at 0 K
        Loss_Options.(Structure).Gap.Value = 0; % Negative value:
        Loss_Options.(Structure).Gap.Weight = 0;
        Loss_Options.(Structure).Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
        Loss_Options.(Structure).Gap.Ref = 'Rocksalt';
    end

end