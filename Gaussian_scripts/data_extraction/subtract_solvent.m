function M_dat_Corrected = subtract_solvent(M_dat,H2O,num_H2O)

    M_dat_Corrected.Pot_E = M_dat.Pot_E - H2O.Pot_E*num_H2O;
    M_dat_Corrected.ZPVE = M_dat.ZPVE - H2O.ZPVE*num_H2O;
    M_dat_Corrected.Thermal_E = M_dat.Thermal_E - H2O.Thermal_E*num_H2O;
    M_dat_Corrected.Thermal_G = M_dat.Thermal_G - H2O.Thermal_G*num_H2O;
    M_dat_Corrected.Total_G = M_dat.Total_G - H2O.Total_G*num_H2O;
    M_dat_Corrected.Trans_E = M_dat.Trans_E - H2O.Trans_E*num_H2O;
    M_dat_Corrected.Trans_TS = M_dat.Trans_TS - H2O.Trans_TS*num_H2O;
    M_dat_Corrected.Rot_E = M_dat.Rot_E - H2O.Rot_E*num_H2O;
    M_dat_Corrected.Rot_TS = M_dat.Rot_TS - H2O.Rot_TS*num_H2O;
    M_dat_Corrected.Vib_E = M_dat.Vib_E - H2O.Vib_E*num_H2O;
    M_dat_Corrected.Vib_TS = M_dat.Vib_TS - H2O.Vib_TS*num_H2O;
    M_dat_Corrected.Tot_E = M_dat.Tot_E - H2O.Tot_E*num_H2O;
    M_dat_Corrected.Tot_TS = M_dat.Tot_TS - H2O.Tot_TS*num_H2O;

end