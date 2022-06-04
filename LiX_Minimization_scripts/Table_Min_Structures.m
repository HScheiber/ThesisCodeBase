Settings = Initialize_MD_Settings;
Data = load(fullfile(Settings.home,'data','MX_JCTF_Min_Data.mat'),'Data').Data;
Theories = {'TF' 'JC' 'JC3P' 'JC4P' 'JCSD'};
%Salts = {'NaCl'};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Structures = {'Rocksalt' 'Wurtzite' 'BetaBeO' 'CsCl' 'FiveFive' 'NiAs' 'AntiNiAs' 'Sphalerite'};

Prop_of_intr = 'E';

Y = nan(length(Theories),length(Salts),length(Structures));
for idx = 1:length(Salts)
    Salt = Salts{idx};
    
    for kdx = 1:length(Theories)
        Theory = Theories{kdx};
        if isfield(Data.(Salt),Theory)
            for jdx = 1:length(Structures)
                Structure = Structures{jdx};
                if isfield(Data.(Salt).(Theory),Structure)
                    Y(kdx,idx,jdx) = Data.(Salt).(Theory).(Structure).(Prop_of_intr);
                end
            end
        end
    end
end


for kdx = 1:length(Theories)
    disp(repmat('*',1,50))
    disp(Theories{kdx})
    disp(repmat('*',1,50))
    
    Y_c = squeeze(Y(kdx,:,:));
    
    disp(num2str(Y_c,'%.2f\t'))
end




