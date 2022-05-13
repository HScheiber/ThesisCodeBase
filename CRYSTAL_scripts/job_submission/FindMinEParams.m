function [a,b,c,Coordinates,run_min] = FindMinEParams(Salt,Structure,Ref_Theory,Theory,Ref_Basis,Basis,...
    home,Ref_Dispersion,Dispersion,Ref_gCP,gCP,Relaxed,DataTypes,P1_Symmetry)

% Initialize the minimization switch as false
run_min = false;

% Load total energies data
X = load([home filesep 'data' filesep 'CRYSTAL' filesep ...
  Salt '_' Structure '_Total_Energies.mat'],'Data');
Data = X.Data;
clearvars X

% Prepare labels
X=1;
Labels{X} = [Salt Label_replace(Structure) '_' Theory Dispersion gCP];
if ~strcmp(Theory,Ref_Theory)
    X=X+1;
    Labels{X} = [Salt Label_replace(Structure) '_' Ref_Theory Dispersion gCP];
    if ~strcmp(Dispersion,Ref_Dispersion)
        X=X+1;
        Labels{X} = [Salt Label_replace(Structure) '_' Ref_Theory Ref_Dispersion gCP];
        if ~strcmp(gCP,Ref_gCP)
            X=X+1;
            Labels{X} = [Salt Label_replace(Structure) '_' Ref_Theory Ref_Dispersion Ref_gCP];
        end
    end
    if ~strcmp(gCP,Ref_gCP)
        X=X+1;
        Labels{X} = [Salt Label_replace(Structure) '_' Ref_Theory Dispersion Ref_gCP];
    end
end
if ~strcmp(Dispersion,Ref_Dispersion)
    X=X+1;
    Labels{X} = [Salt Label_replace(Structure) '_' Theory Ref_Dispersion gCP];
    if ~strcmp(gCP,Ref_gCP)
        X=X+1;
        Labels{X} = [Salt Label_replace(Structure) '_' Theory Ref_Dispersion Ref_gCP];
    end
end
if ~strcmp(gCP,Ref_gCP)
    X=X+1;
    Labels{X} = [Salt Label_replace(Structure) '_' Theory Dispersion Ref_gCP];
end

for j = 1:X
    
    Label = Labels{j};    
    
    % Data of interest if available
    if ~isfield(Data,Label)
        continue
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
    if isempty(DOI_dat) % Try relaxed basis set
        for i = 1:size(DOI,1)
            if strcmp(DOI{i,1},Ref_Basis)
                DOI_dat = DOI{i,2};
                break
            end
        end
    end

    if isempty(DOI_dat)
        continue
    end

    % Find only chosen data type(s) unless relaxed switch active
    Conv_Data = cell(1,9);
    Conv_Ind = ismember(DOI_dat{1,10},DataTypes);
    for i = 2:10
        Conv_Data{i-1} = DOI_dat{1,i}(Conv_Ind);
    end

    % If no converged data, try relaxing constraint 
    if isempty(Conv_Data{1,1}) && Relaxed
        Conv_Data = cell(1,9);
        Conv_Ind = ~ismember(DOI_dat{1,10},DataTypes);
        for i = 2:10
            Conv_Data{i-1} = DOI_dat{1,i}(Conv_Ind);
        end
    end

    % Still no data found, continue
    if isempty(Conv_Data{1,1})
        continue
    end

    % Find lowest energy
    [~,ind] = min(Conv_Data{1,7});

    a = Conv_Data{1,1}(ind);
    b = Conv_Data{1,2}(ind);
    c = Conv_Data{1,3}(ind);
    CoordsCell = Conv_Data{1,6}{ind};

    % Format coordinates
    fcnwrapper = @(x) num2str(x,'%6.12E');
    CoordsCell(:,2:4) = cellfun(fcnwrapper,CoordsCell(:,2:4),'UniformOutput',0);

    % Convert format to crystal input
    for jdx = 1:2
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'LI','3');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'NA','11');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'K','19');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'RB','237');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'CS','255');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'F','9');
        CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'CL','17');
        if strcmpi(Basis,'7-311G')
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'BR','235');
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'I','253');
        elseif strcmpi(Basis,'STO-3G')
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'BR','35');
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'I','53');
        else
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'BR','35');
            CoordsCell{jdx,1} = strrep(CoordsCell{jdx,1},'I','253');
        end
    end
    Coordinates = ['2' newline ...
        CoordsCell{1,1} ' ' CoordsCell{1,2} ' ' CoordsCell{1,3} ' ' CoordsCell{1,4} newline ...
        CoordsCell{2,1} ' ' CoordsCell{2,2} ' ' CoordsCell{2,3} ' ' CoordsCell{2,4}];

    if P1_Symmetry
        MFC = str2double(CoordsCell(1,2:end));
        HFC = str2double(CoordsCell(2,2:end));
        Crystal_Info = P1_Setting(Structure,MFC,HFC);
        Coordinates = num2str(Crystal_Info.N);
        for i = 1:Crystal_Info.N/2
            Coordinates = [Coordinates newline CoordsCell{1,1} ' ' regexprep(num2str(Crystal_Info.FC_Metal(i,:),12),' +',' ')];
            Coordinates = [Coordinates newline CoordsCell{2,1} ' ' regexprep(num2str(Crystal_Info.FC_Halide(i,:),12),' +',' ')];
        end
    end
    return
end

% If nothing found
disp(['No fully minimized data for ' Label '. Running minimization from given inputs.'])
a = [];
b = [];
c = [];
Coordinates = [];
run_min = true;
return

end
