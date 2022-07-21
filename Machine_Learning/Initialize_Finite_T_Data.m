function Finite_T_Data = Initialize_Finite_T_Data(Settings)

    switch Settings.Salt
        case {'CsBr' 'CsI'}
            Finite_T_Data.Structure = 'CsCl';
        otherwise
            Finite_T_Data.Structure = 'Rocksalt';
    end
    
    Finite_T_Data.Fusion_dH = nan;
    Finite_T_Data.Fusion_dV = nan;
    Finite_T_Data.Liquid_V_MP = nan;
    Finite_T_Data.Solid_V_MP = nan;
    Finite_T_Data.MP = nan;
    
    % Load experimental data
    Exp = Load_Experimental_Data;
    Finite_T_Data.Exp_Fusion_dH = Exp.(Settings.Salt).dH_fus; % kJ/mol of Formula units
    Finite_T_Data.Exp_Fusion_dV = Exp.(Settings.Salt).Liquid.V_mp ...
        - Exp.(Settings.Salt).(Finite_T_Data.Structure).V_mp; % Volume change on melting at MP: Ang^3 / Forumla unit
    Finite_T_Data.Exp_Liquid_V_MP = Exp.(Settings.Salt).Liquid.V_mp; % Liquid Volume at MP in Ang^3 / Forumla unit
    Finite_T_Data.Exp_Solid_V_MP = Exp.(Settings.Salt).(Finite_T_Data.Structure).V_mp; % Solid Volume at MP in Ang^3 / Forumla unit
    Finite_T_Data.Exp_MP = Exp.(Settings.Salt).mp; % MP in K
    Finite_T_Data.Exp_Solid_V0 = Exp.(Settings.Salt).(Finite_T_Data.Structure).V_zero;
    
end