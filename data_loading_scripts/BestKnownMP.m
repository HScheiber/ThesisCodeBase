function [T_MP,Matchfound] = BestKnownMP(Settings)

Lit = Load_Literature_Model_MPs;
datasources = fields(Lit);
T_MP = nan;
dmp = Inf;
Matchfound = false;
for idx = 1:length(datasources)
    datasource = datasources{idx};
    if isfield(Lit.(datasource),Settings.Theory)
        if isfield(Lit.(datasource).(Settings.Theory),Settings.Salt)
            if isfield(Lit.(datasource).(Settings.Theory).(Settings.Salt),Settings.Structure)
                if dmp > Lit.(datasource).(Settings.Theory).(Settings.Salt).(Settings.Structure).dmp
                    T_MP = Lit.(datasource).(Settings.Theory).(Settings.Salt).(Settings.Structure).mp;
                    dmp = Lit.(datasource).(Settings.Theory).(Settings.Salt).(Settings.Structure).dmp;
                    Matchfound = true;
                end
            end
        end
    end
end

if isnan(T_MP)
    Exp = Load_Experimental_Data;
    T_MP = Exp.(Settings.Salt).mp;
end


end