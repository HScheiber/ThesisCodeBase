function Structure_out = Find_Ref_Params(Salt,Structure,Theory,Basis,DataType)

    % Load total energies data
    ref_file = [find_home filesep 'data' filesep 'CRYSTAL' filesep ...
      Salt '_' Structure '_Total_Energies.mat'];
  
    if ~isfile(ref_file)
        Structure_out = initialize_crystal(Salt,Structure);
        return
    end
  
    X = load(ref_file,'Data');
    Data = X.Data;
    clearvars X

    % Prepare label
    Label = [Salt Label_replace(Structure) '_' strrep(Theory,'-','_')];

    if ~isfield(Data,Label)
        Structure_out = initialize_crystal(Salt,Structure);
        return
    else
        DOI = Data.(Label);
    end

    % Find basis set match. Try true basis set first, then relax
    DOI_dat = [];
    for i = 1:size(DOI,1)
        if strcmp(DOI{i,1},Basis)
            DOI_dat = DOI{i,2};
            break
        end
    end
    if isempty(DOI_dat)
        Structure_out = initialize_crystal(Salt,Structure);
        return
    end

    % Find only chosen data type(s) unless relaxed switch active
    Conv_Data = cell(1,9);
    Conv_Ind = ismember(DOI_dat{1,10},DataType);
    for i = 2:10
        Conv_Data{i-1} = DOI_dat{1,i}(Conv_Ind);
    end

    % Still no data found, continue
    if isempty(Conv_Data{1,1})
        Structure_out = initialize_crystal(Salt,Structure);
        return
    end

    % Find lowest energy
    [~,ind] = min(Conv_Data{1,7});

    a = Conv_Data{1,1}(ind);
    b = Conv_Data{1,2}(ind);
    c = Conv_Data{1,3}(ind);
    BL = Conv_Data{1,4}(ind);

    % Scale to match the CP2K Input structures (primitive unit cells)
    switch lower(Structure)
        case 'rocksalt'
            Structure_out.a = a/sqrt(2);
            Structure_out.b = b/sqrt(2);
            Structure_out.c = c/sqrt(2);
            Structure_out.SF = a/sqrt(2);
        case 'wurtzite'
            Structure_out.a = a;
            Structure_out.b = b;
            Structure_out.c = c;
            Structure_out.SF = a;
        case 'fivefive'
            Structure_out.a = BL*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Structure_out.b = BL*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
            Structure_out.c = 2*BL;
            Structure_out.SF = BL*1/(sqrt((1/4) + (1/9)*(sind(60)^2)));
        case 'cscl'
            Structure_out.a = a;
            Structure_out.b = b;
            Structure_out.c = c;
            Structure_out.SF = a;
        case 'betabeo'
            Structure_out.a = a;
            Structure_out.b = b;
            Structure_out.c = c;
            Structure_out.SF = a;
        case 'sphalerite'
            Structure_out.a = a/sqrt(2);
            Structure_out.b = b/sqrt(2);
            Structure_out.c = c/sqrt(2);
            Structure_out.SF = a/sqrt(2);
        case 'nias'
            Structure_out.a = a;
            Structure_out.b = b;
            Structure_out.c = c;
            Structure_out.SF = a;
        case 'antinias'
            Structure_out.a = a;
            Structure_out.b = b;
            Structure_out.c = c;
            Structure_out.SF = a;
        otherwise
            error(['Unknown structure: ' Structure])
    end
end
