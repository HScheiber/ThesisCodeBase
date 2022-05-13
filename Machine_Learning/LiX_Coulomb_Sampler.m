function Loss = LiX_Coulomb_Sampler(Charge)
    Salts = {'LiF'};
    Structures = {'Rocksalt' 'NiAs' 'Wurtzite'};
    Conventional = false;
    Model = 'HS';
    Scale = Init_Scaling_Object;
    Scale.r_d.All = 0;
    Scale.b.All = 0;
    Scale.Q = Charge;
    C6_Damp = Init_C6Damping_Object;
    CR_Damp = Init_CRDamping_Object;
    Parallel_gmx = false;
    [GAdjust_MX,GAdjust_MM,GAdjust_XX] = Init_GAdjust_Object;
    
    DFT = Load_Best_DFT_Data('Conventional',Conventional);
    E_Conv = 2625.4996394799; % Hartree -> kJ/mol conversion factor
    
    Loss = 0;
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        for jdx = 1:length(Structures)
            Structure = Structures{jdx};
            Settings = Initialize_MD_Settings;
            Settings.Salt = Salt;
            Settings.Structure = Structure;
            Settings.Use_Conv_cell = Conventional;
            Crystal = Default_Crystal(Settings);
            Crystal.a = DFT.(Salt).(Structure).a;
            Crystal.b = DFT.(Salt).(Structure).b;
            Crystal.c = DFT.(Salt).(Structure).c;
            
            
            Coulomb.(Salt).(Structure) = Sample_Model_energy(Crystal,Salt,Structure,Model,...
                C6_Damp,CR_Damp,Scale,Parallel_gmx,GAdjust_MX,GAdjust_MM,GAdjust_XX,...
                'verbose',false);
            disp([Salt ' ' Structure ' Coulomb Energy: ' num2str(Coulomb.(Salt).(Structure).Coulomb)])
        end
        
        E_Diff_DFT = (DFT.(Salt).(Structures{1}).Non_Interacting_E - ...
            DFT.(Salt).(Structures{2}).Non_Interacting_E)*E_Conv;
        
        E_Diff_Coulomb = (Coulomb.(Salt).(Structures{1}).Coulomb - Coulomb.(Salt).(Structures{2}).Coulomb);
        Loss = Loss + ((E_Diff_DFT - E_Diff_Coulomb)^2);
    end
end
