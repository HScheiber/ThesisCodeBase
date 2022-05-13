function LoadbarObj = PrintCSV(LoadbarObj,Struct_Data,Ymin,filename,Basis_Set)

    N_Col = size(Struct_Data,1);
    Species = Struct_Data(:,4);
    
    Structure_types = Struct_Data(:,2);
    Theory_types = strrep(regexprep(Struct_Data(:,1),'_','-'),'-D0P00000',' (No Dispersion)');
    Energy_min = cell2mat(Ymin(:));
    a_min = cell2mat(Struct_Data(:,3));
    b_min = cell2mat(Struct_Data(:,5));
    c_min = cell2mat(Struct_Data(:,6));
    aT_min = cell2mat(Struct_Data(:,14));
    bT_min = cell2mat(Struct_Data(:,15));
    cT_min = cell2mat(Struct_Data(:,16));
    
    FC_Mx = nan(N_Col,1);
    FC_My = nan(N_Col,1);
    FC_Mz = nan(N_Col,1);
    FC_Xx = nan(N_Col,1);
    FC_Xy = nan(N_Col,1);
    FC_Xz = nan(N_Col,1);
    
    
    Metal_Basis_Sets = cell(N_Col,1);
    Halide_Basis_Sets = cell(N_Col,1);
    
    Empirical_models = ~cellfun(@isempty,regexp(Theory_types,'(JC)|(TF)'));
    
    switch Basis_Set{1}
        case 'Mpob-TZVP_HECP28MDF-VTZ'
            Metal_Basis_Sets(:) = {'pob-TZVP'};
            Halide_Basis_Sets(:) = {'ECP28MDF-VTZ'};
        otherwise
            Metal_Basis_Sets(:) = Basis_Set;
            Halide_Basis_Sets(:) = Basis_Set;
    end
    Metal_Basis_Sets(Empirical_models) = {''};
    Halide_Basis_Sets(Empirical_models) = {''};
    
    for i=1:N_Col
        FC = Struct_Data{i,9};
        
        if any(strcmpi(FC{1,1},{'Li' 'Na' 'K' 'Rb' 'Cs'}))
            FC_Metal = [FC{1,2:end}];
            FC_Halide = [FC{2,2:end}];
        else
            FC_Metal = [FC{2,2:end}];
            FC_Halide = [FC{1,2:end}];
        end
        
        
        FC_Mx(i) = mod(FC_Metal(1),1);
        FC_My(i) = mod(FC_Metal(2),1);
        FC_Mz(i) = mod(FC_Metal(3),1);
        FC_Xx(i) = mod(FC_Halide(1),1);
        FC_Xy(i) = mod(FC_Halide(2),1);
        FC_Xz(i) = mod(FC_Halide(3),1);
    end

    
    waitbar(0.25,LoadbarObj,'Generating Table...')

    % Merge into table
    Data_Table = table(Species,Structure_types,Theory_types,Metal_Basis_Sets,Halide_Basis_Sets,...
        Energy_min,a_min,b_min,c_min,aT_min,bT_min,cT_min,FC_Mx,FC_My,FC_Mz,FC_Xx,FC_Xy,FC_Xz,...
        'VariableNames',{'Species' 'Structure' 'Theory' 'Metal_Basis' 'Halide_Basis'...
        'E_min_kJmol' 'a_Ang' 'b_Ang' 'c_Ang' 'aT_Ang' 'bT_Ang' 'cT_Ang' 'Mx' 'My' 'Mz' 'Xx' 'Xy' 'Xz'});
    
    % Sort Table
    Data_Table = sortrows(Data_Table,{'Species' 'Theory' 'Structure'},{'ascend','ascend','ascend'});   
    writetable(Data_Table,filename)
 
end