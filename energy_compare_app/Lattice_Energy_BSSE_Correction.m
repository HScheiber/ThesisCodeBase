% Function to correct for Basis set superposition error
function Data = Lattice_Energy_BSSE_Correction(Data,Metal_BSSE_Data,Halide_BSSE_Data,Label)

ST = regexp(Label,'(?<=.{2,4}_).','match');
Salt_Type = ST{1};
tolerance = 0.01; % Tolerance for numerical comparisons

Fields = fieldnames(Data);
M = length(Fields);
for i = 1:M % Loop through calculation types
    
    % Current calculation type
    calc_type = regexprep(strrep(Fields{i},[Label '_'],''),'_D[2-3].*','');
    
    % Number of basis sets
    N = size(Data.(Fields{i}),1); 
    for j = N:-1:1 % Loop through basis sets
        matched_Metal = false;
        matched_Halide = false;
        
        % Current Basis set
        Basis_set = Data.(Fields{i}){j,1};
        
        % Remove if data is missing
        if isempty(Data.(Fields{i}){j,2}{3}) || isempty(Data.(Fields{i}){j,2}{1,2})
            Data.(Fields{i})(j,:) = [];
            continue
        end
        
        % Correct for Metal first
        if ~isempty(Metal_BSSE_Data)
            Metal_Fields = fieldnames(Metal_BSSE_Data);
            O = length(Metal_Fields);
            for k = 1:O

                Metal_calc_type = strrep(Metal_Fields{k},[Label '_'],'');

                % Check if salt type and calculation type matches
                if isempty(strfind(Metal_Fields{k},Label)) || ~strcmp(calc_type,Metal_calc_type)
                    continue
                end

                % Number of Metal basis sets
                P = size(Metal_BSSE_Data.(Metal_Fields{k}),1);
                for m = 1:P

                    % Metal basis set
                    Metal_Basis_set = Metal_BSSE_Data.(Metal_Fields{k}){m,1};

                    % Check for basis set match (pob-TZVP with def2-TZVPD
                    % considered a match)
                    if strcmp(Basis_set,'pob-TZVP') && strcmp(Metal_Basis_set,'def2-TZVPD')
                        matched_Metal = true;
                    elseif strcmp(Basis_set,'pob-TZVP') && strcmp(Metal_Basis_set,'pob-TZVP')
                        continue
                    elseif strcmp(Basis_set,'def2-TZVPD') && strcmp(Metal_Basis_set,'def2-TZVPD')
                        continue
                    elseif strcmp(Basis_set,Metal_Basis_set)
                        matched_Metal = true;
                    else
                        continue
                    end

                    % Loop through PES
                    Q = length(Data.(Fields{i}){j,2}{1,2});
                    R = length(Metal_BSSE_Data.(Metal_Fields{k}){m,2}{1,2});
                    for n = Q:-1:1
                        matched = false;
                        for o = 1:R
                            if abs(Data.(Fields{i}){j,2}{1,2}(n) - Metal_BSSE_Data.(Metal_Fields{k}){m,2}{1,2}(o)) < tolerance
                                Data.(Fields{i}){j,2}{1,8}(n) = Data.(Fields{i}){j,2}{1,8}(n) - Metal_BSSE_Data.(Metal_Fields{k}){m,2}{1,8}(o);
                                matched = true;
                                break
                            end
                        end
                        if ~matched
                            point = num2str(Data.(Fields{i}){j,2}{1,2}(n));
                            Data.(Fields{i}){j,2}{1,2}(n) = [];
                            Data.(Fields{i}){j,2}{1,3}(n) = [];
                            Data.(Fields{i}){j,2}{1,4}(n) = [];
                            Data.(Fields{i}){j,2}{1,5}(n) = [];
                            Data.(Fields{i}){j,2}{1,6}(n) = [];
                            Data.(Fields{i}){j,2}{1,7}(n) = [];
                            Data.(Fields{i}){j,2}{1,8}(n) = [];
                            Data.(Fields{i}){j,2}{1,9}(n) = [];
                            Data.(Fields{i}){j,2}{1,10}(n) = [];
                        end
                    end
                    if matched_Metal
                        break
                    end
                end
                if matched_Metal
                    break
                end
            end
        end
        
        % Correct for Halide second
        if ~isempty(Halide_BSSE_Data)
            Halide_Fields = fieldnames(Halide_BSSE_Data);
            O = length(Halide_Fields);
            for k = 1:O

                Halide_calc_type = strrep(Halide_Fields{k},[Label '_'],'');

                % Check if salt type matches
                if isempty(strfind(Halide_Fields{k},Label)) || ~strcmp(calc_type,Halide_calc_type)
                    continue
                end

                % Number of Halide basis sets
                P = size(Halide_BSSE_Data.(Halide_Fields{k}),1);
                for m = 1:P

                    % Halide basis set
                    Halide_Basis_set = Halide_BSSE_Data.(Halide_Fields{k}){m,1};

                    % Check for basis set match (pob-TZVP with def2-TZVPD
                    % considered a match)
                    if strcmp(Basis_set,'pob-TZVP') && strcmp(Halide_Basis_set,'def2-TZVPD')
                        matched_Halide = true;
                    elseif strcmp(Basis_set,'pob-TZVP') && strcmp(Halide_Basis_set,'pob-TZVP')
                        continue
                    elseif strcmp(Basis_set,'def2-TZVPD') && strcmp(Halide_Basis_set,'def2-TZVPD')
                        continue
                    elseif strcmp(Basis_set,Halide_Basis_set)
                        matched_Halide = true;
                    else
                        continue
                    end

                    % Loop through PES
                    Q = length(Data.(Fields{i}){j,2}{1,2});
                    R = length(Halide_BSSE_Data.(Halide_Fields{k}){m,2}{1,2});
                    for n = Q:-1:1
                        matched = false;
                        for o = 1:R
                            if abs(Data.(Fields{i}){j,2}{1,2}(n) - Halide_BSSE_Data.(Halide_Fields{k}){m,2}{1,2}(o)) < tolerance
                                Data.(Fields{i}){j,2}{1,8}(n) = Data.(Fields{i}){j,2}{1,8}(n) - ...
                                    Halide_BSSE_Data.(Halide_Fields{k}){m,2}{1,8}(o);
                                matched = true;
                                break
                            end
                        end
                        if ~matched
                            point = num2str(Data.(Fields{i}){j,2}{1,2}(n));
                            Data.(Fields{i}){j,2}{1,2}(n) = [];
                            Data.(Fields{i}){j,2}{1,3}(n) = [];
                            Data.(Fields{i}){j,2}{1,4}(n) = [];
                            Data.(Fields{i}){j,2}{1,5}(n) = [];
                            Data.(Fields{i}){j,2}{1,6}(n) = [];
                            Data.(Fields{i}){j,2}{1,7}(n) = [];
                            Data.(Fields{i}){j,2}{1,8}(n) = [];
                            Data.(Fields{i}){j,2}{1,9}(n) = [];
                            Data.(Fields{i}){j,2}{1,10}(n) = [];
                        end
                    end
                    if matched_Halide
                        break
                    end
                end
                if matched_Halide
                    break
                end
            end
        end
        
        % Prune missing data
        if ~matched_Metal || ~matched_Halide || isempty(Data.(Fields{i}){j,2}{1,2})
            Data.(Fields{i})(j,:) = [];
        end
    end
    % Prune empty fields
    if isempty(Data.(Fields{i}))
        Data = rmfield(Data,Fields{i});
    end
end

end