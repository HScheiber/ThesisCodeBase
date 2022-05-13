%% Inputs
Theories = {'JC'}; %{'HSEsol' 'PW1PW_D3' 'PW1PW_D4' 'JC' 'TF'};
Salts = {'LiCl'}; % 'LiF' 'LiCl' 'LiBr' 'LiI' 
Ref_Structure = 'Wurtzite';
Structures = {'FiveFive' 'Rocksalt'}; %'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'}; % 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'
Merge_Tol = 0.01; % Angstroms
Shell_Tol = 50.0; % Angstroms
if ispc
    datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\';
else
    datadir = '/home/user/ThesisCodeBase/data/';
end
Data_Type = 1;
Basis_Set = 'pob-TZVP';



for ii = 1:length(Salts)
    Salt = Salts{ii};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    
    for jj = 1:length(Theories)
        Theory = Theories{jj};

        
        %% Get reference structure's data
        
        % Find equilibrium lattice parameters of given Salt/Structure/Theory
        if contains(Theory,{'TF' 'JC'}) % empirical models
            load(fullfile(datadir,'GROMACS',[Salt '_' Ref_Structure '_Lattice_Energies.mat']),'Data');

            try
                SubData = Data.(Salt).(Ref_Structure).(Theory);
            catch
                warning(['Missing Data for: ' Salt ' ' Ref_Structure ' ' Theory]);
                continue
            end
            DT = [SubData{:,9}];
            DatInd = ismember(DT,Data_Type);

            if sum(DatInd) == 0
                continue
            end

            CryStruc.a = [SubData{DatInd,1}];
            CryStruc.b = [SubData{DatInd,2}];
            CryStruc.c = [SubData{DatInd,3}];
            CellAngles = SubData(DatInd,8);
            CryStruc.FC = SubData(DatInd,6);   
            if sum(DatInd) > 1
                En = [SubData{DatInd,7}];
                [~,midx] = min(En);
                CryStruc.a = CryStruc.a(midx);
                CryStruc.b = CryStruc.b(midx);
                CryStruc.c = CryStruc.c(midx);
                CryStruc.alpha = CellAngles{midx}(1);
                CryStruc.beta = CellAngles{midx}(2);
                CryStruc.gamma = CellAngles{midx}(3);
                CryStruc.FC = CryStruc.FC{midx};
            else
                CryStruc.alpha = CellAngles{1}(1);
                CryStruc.beta = CellAngles{1}(2);
                CryStruc.gamma = CellAngles{1}(3);
                CryStruc.FC = CryStruc.FC{1};
            end
        else % DFT data
            load(fullfile(datadir,'CRYSTAL',[Salt '_' Ref_Structure '_Lattice_Energies.mat']),'Data');

            Label = [Salt '_' Ref_Structure(1) '_' Theory];
            SubData = Data.(Label);

            DOI = [];
            for i = 1:size(SubData,1)
                if strcmp(SubData{i,1},Basis_Set)
                    DOI = SubData{i,2};
                    break
                end
            end
            if ~isempty(DOI)
                DT = DOI{10};

                DatInd = ismember(DT,Data_Type);

                if sum(DatInd) == 0
                    continue
                end

                CryStruc.a = DOI{2}(DatInd);
                CryStruc.b = DOI{3}(DatInd);
                CryStruc.c = DOI{4}(DatInd);
                CryStruc.alpha = DOI{11}(DatInd);
                CryStruc.beta = DOI{12}(DatInd);
                CryStruc.gamma = DOI{13}(DatInd);
                CryStruc.FC = DOI{7}(DatInd);        
                if sum(DatInd) > 1
                    En = DOI{8}(DatInd);
                    [~,midx] = min(En);
                    CryStruc.a = CryStruc.a(midx);
                    CryStruc.b = CryStruc.b(midx);
                    CryStruc.c = CryStruc.c(midx);
                    CryStruc.alpha = CryStruc.alpha(midx);
                    CryStruc.beta = CryStruc.beta(midx);
                    CryStruc.gamma = CryStruc.gamma(midx);
                    CryStruc.FC = CryStruc.FC{midx};
                else
                    CryStruc.FC = CryStruc.FC{1};
                end
            else
                continue
            end
        end

        % Generate coordination shells using structural parameters
        [MX_BL,MM_BL,XX_BL,XM_BL] = Gen_Stars(CryStruc,Ref_Structure,Metal,Halide,10);

        % Collapse coordination shells if less than the specified merge tolerance
        for ix = size(MX_BL,1):-1:2
            if abs(MX_BL(ix,1) - MX_BL(ix-1,1)) < Merge_Tol
                MX_BL(ix-1,1) = (MX_BL(ix,1)*MX_BL(ix,2) + MX_BL(ix-1,1)*MX_BL(ix-1,2))/...
                    (MX_BL(ix,2)+MX_BL(ix-1,2)); % new distance = weighted mean
                MX_BL(ix-1,2) = MX_BL(ix-1,2) + MX_BL(ix,2); % New CN is sum of old
                MX_BL(ix,:) = [];
            end
        end
        for ix = size(MM_BL,1):-1:2
            if abs(MM_BL(ix,1) - MM_BL(ix-1,1)) < Merge_Tol
                MM_BL(ix-1,1) = (MM_BL(ix,1)*MM_BL(ix,2) + MM_BL(ix-1,1)*MM_BL(ix-1,2))/...
                    (MM_BL(ix,2)+MM_BL(ix-1,2)); % new distance = weighted mean
                MM_BL(ix-1,2) = MM_BL(ix-1,2) + MM_BL(ix,2); % New CN is sum of old
                MM_BL(ix,:) = [];
            end
        end
        for ix = size(XX_BL,1):-1:2
            if abs(XX_BL(ix,1) - XX_BL(ix-1,1)) < Merge_Tol
                XX_BL(ix-1,1) = (XX_BL(ix,1)*XX_BL(ix,2) + XX_BL(ix-1,1)*XX_BL(ix-1,2))/...
                    (XX_BL(ix,2)+XX_BL(ix-1,2)); % new distance = weighted mean
                XX_BL(ix-1,2) = XX_BL(ix-1,2) + XX_BL(ix,2); % New CN is sum of old
                XX_BL(ix,:) = [];
            end
        end

        MX_BL_A = MX_BL(MX_BL(1,:) <= MX_BL(1,1)+Shell_Tol,1);
        MX_CN_A = MX_BL(MX_BL(1,:) <= MX_BL(1,1)+Shell_Tol,2);
        XX_BL_A = XX_BL(XX_BL(1,:) <= XX_BL(1,1)+Shell_Tol,1);
        XX_CN_A = XX_BL(XX_BL(1,:) <= XX_BL(1,1)+Shell_Tol,2);
        MM_BL_A = MM_BL(MM_BL(1,:) <= MM_BL(1,1)+Shell_Tol,1);
        MM_CN_A = MM_BL(MM_BL(1,:) <= MM_BL(1,1)+Shell_Tol,2);
        
        C6_MX_CA = sum(MX_CN_A./(MX_BL_A.^6));
        C6_XX_CA = sum(XX_CN_A./(XX_BL_A.^6));
        C6_MM_CA = sum(MM_CN_A./(MM_BL_A.^6));
        C8_MX_CA = sum(MX_CN_A./(MX_BL_A.^8));
        C8_XX_CA = sum(XX_CN_A./(XX_BL_A.^8));
        C8_MM_CA = sum(MM_CN_A./(MM_BL_A.^8));

        %% Get comparison structure's data
        for iii = 1:length(Structures)
            Structure = Structures{iii};

            % Find equilibrium lattice parameters of given Salt/Structure/Theory
            if contains(Theory,{'TF' 'JC'}) % empirical models
                load(fullfile(datadir,'GROMACS',[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
                
                try
                    SubData = Data.(Salt).(Structure).(Theory);
                catch
                    warning(['Missing Data for: ' Salt ' ' Structure ' ' Theory]);
                    continue
                end
                DT = [SubData{:,9}];
                DatInd = ismember(DT,Data_Type);

                if sum(DatInd) == 0
                    continue
                end

                CryStruc.a = [SubData{DatInd,1}];
                CryStruc.b = [SubData{DatInd,2}];
                CryStruc.c = [SubData{DatInd,3}];
                CellAngles = SubData(DatInd,8);
                CryStruc.FC = SubData(DatInd,6);   
                if sum(DatInd) > 1
                    En = [SubData{DatInd,7}];
                    [~,midx] = min(En);
                    CryStruc.a = CryStruc.a(midx);
                    CryStruc.b = CryStruc.b(midx);
                    CryStruc.c = CryStruc.c(midx);
                    CryStruc.alpha = CellAngles{midx}(1);
                    CryStruc.beta = CellAngles{midx}(2);
                    CryStruc.gamma = CellAngles{midx}(3);
                    CryStruc.FC = CryStruc.FC{midx};
                else
                    CryStruc.alpha = CellAngles{1}(1);
                    CryStruc.beta = CellAngles{1}(2);
                    CryStruc.gamma = CellAngles{1}(3);
                    CryStruc.FC = CryStruc.FC{1};
                end
            else % DFT data
                load(fullfile(datadir,'CRYSTAL',[Salt '_' Structure '_Lattice_Energies.mat']),'Data');

                Label = [Salt '_' Structure(1) '_' Theory];
                SubData = Data.(Label);

                DOI = [];
                for i = 1:size(SubData,1)
                    if strcmp(SubData{i,1},Basis_Set)
                        DOI = SubData{i,2};
                        break
                    end
                end
                if ~isempty(DOI)
                    DT = DOI{10};

                    DatInd = ismember(DT,Data_Type);

                    if sum(DatInd) == 0
                        continue
                    end

                    CryStruc.a = DOI{2}(DatInd);
                    CryStruc.b = DOI{3}(DatInd);
                    CryStruc.c = DOI{4}(DatInd);
                    CryStruc.alpha = DOI{11}(DatInd);
                    CryStruc.beta = DOI{12}(DatInd);
                    CryStruc.gamma = DOI{13}(DatInd);
                    CryStruc.FC = DOI{7}(DatInd);        
                    if sum(DatInd) > 1
                        En = DOI{8}(DatInd);
                        [~,midx] = min(En);
                        CryStruc.a = CryStruc.a(midx);
                        CryStruc.b = CryStruc.b(midx);
                        CryStruc.c = CryStruc.c(midx);
                        CryStruc.alpha = CryStruc.alpha(midx);
                        CryStruc.beta = CryStruc.beta(midx);
                        CryStruc.gamma = CryStruc.gamma(midx);
                        CryStruc.FC = CryStruc.FC{midx};
                    else
                        CryStruc.FC = CryStruc.FC{1};
                    end
                else
                    continue
                end
            end

            % Generate coordination shells using structural parameters
            [MX_BL,MM_BL,XX_BL,XM_BL] = Gen_Stars(CryStruc,Structure,Metal,Halide,10);

            % Collapse coordination shells if less than the specified merge tolerance
            for ix = size(MX_BL,1):-1:2
                if abs(MX_BL(ix,1) - MX_BL(ix-1,1)) < Merge_Tol
                    MX_BL(ix-1,1) = (MX_BL(ix,1)*MX_BL(ix,2) + MX_BL(ix-1,1)*MX_BL(ix-1,2))/...
                        (MX_BL(ix,2)+MX_BL(ix-1,2)); % new distance = weighted mean
                    MX_BL(ix-1,2) = MX_BL(ix-1,2) + MX_BL(ix,2); % New CN is sum of old
                    MX_BL(ix,:) = [];
                end
            end
            for ix = size(MM_BL,1):-1:2
                if abs(MM_BL(ix,1) - MM_BL(ix-1,1)) < Merge_Tol
                    MM_BL(ix-1,1) = (MM_BL(ix,1)*MM_BL(ix,2) + MM_BL(ix-1,1)*MM_BL(ix-1,2))/...
                        (MM_BL(ix,2)+MM_BL(ix-1,2)); % new distance = weighted mean
                    MM_BL(ix-1,2) = MM_BL(ix-1,2) + MM_BL(ix,2); % New CN is sum of old
                    MM_BL(ix,:) = [];
                end
            end
            for ix = size(XX_BL,1):-1:2
                if abs(XX_BL(ix,1) - XX_BL(ix-1,1)) < Merge_Tol
                    XX_BL(ix-1,1) = (XX_BL(ix,1)*XX_BL(ix,2) + XX_BL(ix-1,1)*XX_BL(ix-1,2))/...
                        (XX_BL(ix,2)+XX_BL(ix-1,2)); % new distance = weighted mean
                    XX_BL(ix-1,2) = XX_BL(ix-1,2) + XX_BL(ix,2); % New CN is sum of old
                    XX_BL(ix,:) = [];
                end
            end
            
            MX_BL_B = MX_BL(MX_BL(1,:) <= MX_BL(1,1)+Shell_Tol,1);
            MX_CN_B = MX_BL(MX_BL(1,:) <= MX_BL(1,1)+Shell_Tol,2);
            XX_BL_B = XX_BL(XX_BL(1,:) <= XX_BL(1,1)+Shell_Tol,1);
            XX_CN_B = XX_BL(XX_BL(1,:) <= XX_BL(1,1)+Shell_Tol,2);
            MM_BL_B = MM_BL(MM_BL(1,:) <= MM_BL(1,1)+Shell_Tol,1);
            MM_CN_B = MM_BL(MM_BL(1,:) <= MM_BL(1,1)+Shell_Tol,2);
            
            C6_MX_CB = sum(MX_CN_B./(MX_BL_B.^6));
            C6_XX_CB = sum(XX_CN_B./(XX_BL_B.^6));
            C6_MM_CB = sum(MM_CN_B./(MM_BL_B.^6));
            C8_MX_CB = sum(MX_CN_B./(MX_BL_B.^8));
            C8_XX_CB = sum(XX_CN_B./(XX_BL_B.^8));
            C8_MM_CB = sum(MM_CN_B./(MM_BL_B.^8));
            
            C6_MX_Critereon = C6_MX_CA/C6_MX_CB;
            C6_XX_Critereon = C6_XX_CA/C6_XX_CB;
            C6_MM_Critereon = C6_MM_CA/C6_MM_CB;
            C8_MX_Critereon = C8_MX_CA/C8_MX_CB;
            C8_XX_Critereon = C8_XX_CA/C8_XX_CB;
            C8_MM_Critereon = C8_MM_CA/C8_MM_CB;

            disp([Salt ' ' Structure ' ' Theory])
            disp(['Criterion6_MX = ' num2str(C6_MX_Critereon,'%2.2f')])
            disp(['Criterion8_MX = ' num2str(C8_MX_Critereon,'%2.2f')])
            disp(['Criterion6_XX = ' num2str(C6_XX_Critereon,'%2.2f')])
            disp(['Criterion8_XX = ' num2str(C8_XX_Critereon,'%2.2f')])
            disp(['Criterion6_MM = ' num2str(C6_MM_Critereon,'%2.2f')])
            disp(['Criterion8_MM = ' num2str(C8_MM_Critereon,'%2.2f')])
        end
    end
end

