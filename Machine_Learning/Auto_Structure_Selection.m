function Structures = Auto_Structure_Selection(Loss_Options)
    tol = sqrt(eps);
    Loss_fields = fields(Loss_Options);
    Structures = {};
    
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
                        if strcmp(subfield,'RLE') && ~any(cellfun(@(x) strcmp(x,'Rocksalt'),Structures))
                            Structures{end+1} = 'Rocksalt';
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
end