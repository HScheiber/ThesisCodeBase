N_Row = height(Data_Table_Cell_Opt);

for idx = 1:N_Row
    Salt = Data_Table_Cell_Opt(idx,:).Species{1};
    Structure = Data_Table_Cell_Opt(idx,:).Structure{1};
    Theory = Data_Table_Cell_Opt(idx,:).Theory{1};
    Metal_Basis = Data_Table_Cell_Opt(idx,:).Metal_Basis{1};
    Halide_Basis = Data_Table_Cell_Opt(idx,:).Halide_Basis{1};
    
    % Select data from full opt cell
    Salt_idx = strcmp(Salt,Data_Table_Full_Opt.Species);
    Structure_idx = strcmp(Structure,Data_Table_Full_Opt.Structure);
    Theory_idx = strcmp(Theory,Data_Table_Full_Opt.Theory);
    Metal_Basis_idx = strcmp(Metal_Basis,Data_Table_Full_Opt.Metal_Basis);
    Halide_Basis_idx = strcmp(Halide_Basis,Data_Table_Full_Opt.Halide_Basis);
    
    full_idx = (Salt_idx & Structure_idx & Theory_idx & Metal_Basis_idx & Halide_Basis_idx);
    sel_data_full_opt = Data_Table_Full_Opt(full_idx,:);
    
    if isempty(sel_data_full_opt)
        continue
    end
    
    % Check if the structure is correct
    switch Structure
        case 'BetaBeO'
            if (sel_data_full_opt.Mx - 0.25) < 0.03 && (sel_data_full_opt.Xx - 0.25) < 0.03
                continue
            else
                Data_Table_Cell_Opt(idx,:) = sel_data_full_opt; %#ok<*SAGROW>
            end
        case 'Wurtzite'
            if (abs(sel_data_full_opt.Mz - sel_data_full_opt.Xz) - 0.5 < 0.03 )
                continue
            else
                Data_Table_Cell_Opt(idx,:) = sel_data_full_opt; %#ok<*SAGROW>
            end
        case 'FiveFive'
            if (sel_data_full_opt.My - 0.25 < 0.03 ) && (sel_data_full_opt.Xy - 0.25) < 0.03
                continue
            else
                Data_Table_Cell_Opt(idx,:) = sel_data_full_opt; %#ok<*SAGROW>
            end
        otherwise
            Data_Table_Cell_Opt(idx,:) = sel_data_full_opt; %#ok<*SAGROW>
    end
end

writetable(Data_Table_Full_Opt,'LiX_Dataset.csv')
    