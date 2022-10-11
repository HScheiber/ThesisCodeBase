% Function which generates default structure parameters
function Model_Scaled = ModelName(Settings)

    % Parameter set for TF model
    if Settings.TF_Paramset ~= 0 && strcmp(Settings.Theory,'TF')
        Model_Scaled = [Settings.Theory num2str(Settings.TF_Paramset)];
    elseif Settings.TF_Paramset ~= 0 && ~strcmp(Settings.Theory,'TF') % Skip JC models
        Model_Scaled = '';
        return
    else
        Model_Scaled = Settings.Theory;
    end
    
    if Settings.C6_Damp.MX == Settings.C6_Damp.MM && Settings.C6_Damp.MX == Settings.C6_Damp.XX
        if Settings.C6_Damp.MX ~= 0
            Model_Scaled = [Model_Scaled 'd' num2str(Settings.C6_Damp.MX)];
        end
    else
        for current_interaction = fieldnames(Settings.C6_Damp)'
            tag = ['d' current_interaction{1}];

            if ~contained_in_cell(Settings.C6_Damp.(current_interaction{1}),{'MM' 'XX' 'MX' })
                continue
            elseif Settings.C6_Damp.(current_interaction{1}) ~= 0
                Model_Scaled = [Model_Scaled tag num2str(Settings.C6_Damp.(current_interaction{1}))];
            end
        end
    end
    
    for current_scale = fieldnames(Settings.S)'
        
        if strcmp(current_scale{1},'Q')

            tag = ['_' current_scale{1}];
            if ~ismembertol(1.0,Settings.S.(current_scale{1}),1e-2)
                mtxt = strrep(num2str(Settings.S.(current_scale{1}),'%10.2f'),'-','N');
                Model_Scaled = [Model_Scaled tag mtxt];
            end
        elseif strcmp(current_scale{1},'TFParamset')
            continue
        elseif strcmp(current_scale{1},'N')
            continue
        elseif strcmp(current_scale{1},'n')
            continue
        else
        
            for current_interaction = fieldnames(Settings.S.(current_scale{1}))'

                if ( strcmp(current_scale{1},'E') || strcmp(current_scale{1},'S') )&& ...
                        ~ismembertol(1.0,Settings.S.(current_scale{1}).(current_interaction{1}),1e-5) && ...
                        strcmp(Settings.Theory,'TF')
                    Model_Scaled = ''; % epsilon scaling not defined for TF model
                    return
                end

                if strcmp(current_scale{1},'A') && ...
                        ~ismembertol(1.0,Settings.S.(current_scale{1}).(current_interaction{1}),1e-5) && ...
                        strcmp(Settings.Theory,'JC')
                    Model_Scaled = ''; % alpha scaling not defined for JC model
                    return
                end

                if strcmp(current_interaction{1},'All')
                    tag = ['_' current_scale{1}];
                else
                    tag = ['_' current_interaction{1} current_scale{1}];
                end

                if ~ismembertol(1.0,Settings.S.(current_scale{1}).(current_interaction{1}),1e-5)
                    mtxt = strrep(num2str(Settings.S.(current_scale{1}).(current_interaction{1}),'%10.5f'),'-','N');
                    Model_Scaled = [Model_Scaled tag mtxt];
                end
            end
        end
    end

    for current_interaction = fieldnames(Settings.CR_Damp)'
        if strcmp(current_interaction{1},'All')
            tag = '_dCR';
        else
            tag = ['_dCR' current_interaction{1}];
        end

        if Settings.CR_Damp.(current_interaction{1}).r_d >= 0 && Settings.CR_Damp.(current_interaction{1}).b >= 0
            mtxt = ['rd' num2str(Settings.CR_Damp.(current_interaction{1}).r_d,'%10.3f') 'b' num2str(Settings.CR_Damp.(current_interaction{1}).b,'%10.3f')];
            Model_Scaled = [Model_Scaled tag mtxt];
        end
    end
    
    % For naming purposes
    G_Adj = '';
    if sum(Settings.GAdjust_MX(:,1)) ~= 0
        G_Adj = [G_Adj '_GMX'];
    end

    if sum(Settings.GAdjust_MM(:,1)) ~= 0
        G_Adj = [G_Adj '_GMM'];
    end

    if sum(Settings.GAdjust_XX(:,1)) ~= 0
        G_Adj = [G_Adj '_GXX'];
    end
    
    Model_Scaled = [Model_Scaled G_Adj];
    
    % Fixed names for some models
    Model_Scaled = strrep(Model_Scaled,'JC_MMD910.00000_XXD0.23000_dCRMMrd0.260b75.000','JC_ModelA');
    
end