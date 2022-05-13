function MP = Load_Model_MP(Salt,Theory,Model)

    home = find_home;
    DataBase = fullfile(home,'data','Melting_Point_Data.mat');

    Data = load(DataBase,'Data');
    Data = Data.Data;
    
    if isempty(Model)
        Full_Model = Theory;
    else
        Full_Model = [Theory '_Model_' Model];
    end
    
    try
        Calcs = fields(Data.(Salt).(Full_Model).Rocksalt);
        
        N_atoms = zeros(length(Calcs),1);
        for idx = 1:length(Calcs)
            Calc_txt = Calcs{idx};
            if Data.(Salt).(Full_Model).Rocksalt.(Calc_txt).IsBubble
                continue
            else
                N_atoms(idx) = Data.(Salt).(Full_Model).Rocksalt.(Calc_txt).N_Total;
            end
        end
        
        [~,midx] = max(N_atoms);
        
        mps_idx = (~Data.(Salt).(Full_Model).Rocksalt.(Calcs{midx}).Freeze_Trace & ...
            ~Data.(Salt).(Full_Model).Rocksalt.(Calcs{midx}).Melt_Trace);
        if sum(mps_idx) == 0
            MP = mean(Data.(Salt).(Full_Model).Rocksalt.(Calcs{midx}).dT);
        else
            MPs = Data.(Salt).(Full_Model).Rocksalt.(Calcs{midx}).T_Trace(mps_idx);
            MP = mean([min(MPs) max(MPs)]);
        end
        
    catch
        disp('No known melting point: using experimental melting point.')
        Data = Load_Experimental_Data;
        MP = Data.(Salt).mp;
    end

end