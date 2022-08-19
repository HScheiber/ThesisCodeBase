function Structures = Auto_Structure_Selection(Model)
    tol = sqrt(eps);
    Loss_Options = Model.Loss_Options;
    Loss_fields = fields(Loss_Options);
    Structures = {};
    
    switch Model.Salt
        case {'CsBr' 'CsI'}
            Finite_T_Data_Structure = 'CsCl';
        otherwise
            Finite_T_Data_Structure = 'Rocksalt';
    end
    
    
    for idx = 1:length(Loss_fields)
        Loss_field  = Loss_fields{idx};
        skip_remain = false;
        if isstruct(Loss_Options.(Loss_field))
            
            subfields = fields(Loss_Options.(Loss_field));
            for jdx = 1:length(subfields)
                subfield = subfields{jdx};
                
                if isstruct(Loss_Options.(Loss_field).(subfield))
                    if Loss_Options.(Loss_field).(subfield).Weight > tol
                        if ~skip_remain
                            Structures{end+1} = Loss_field;
                        end
                        
                        % Make sure Reference is also present
                        if ~any(cellfun(@(x) strcmp(x,Loss_Options.(Loss_field).(subfield).Ref),Structures))
                            Structures{end+1} = Loss_Options.(Loss_field).(subfield).Ref;
                        end
                        break
                    else
                        continue
                    end
                elseif ~skip_remain
                    % Check if 0
                    if abs(Loss_Options.(Loss_field).(subfield)) > tol
                        Structures{end+1} = Loss_field;
                        
                        % If RLE, make sure reference is also present
                        if strcmp(subfield,'RLE') && ~any(cellfun(@(x) strcmp(x,Finite_T_Data_Structure),Structures))
                            Structures{end+1} = Finite_T_Data_Structure;
                        end
                        skip_remain = true;
                        continue
                    else
                        continue
                    end
                end
            end
        end
    end
    for idx = 1:length(Loss_fields)
        Loss_field  = Loss_fields{idx};
        if contained_in_cell(Loss_field,{'Fusion_Enthalpy', 'MP_Volume_Change', 'Liquid_MP_Volume', 'Solid_MP_Volume', 'MP', 'Liquid_DM_MP'})
            if Loss_Options.(Loss_field) > tol && ~any(cellfun(@(x) strcmp(x,'Rocksalt'),Structures))
                Structures{end+1} = 'Rocksalt';
                break
            end
        end
    end
end