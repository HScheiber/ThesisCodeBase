% Function which holds pre-evaluated Gaussian adjustment parameters
function [GAdjust_MX,GAdjust_MM,GAdjust_XX] = GAjustPreset(JobID,GAdjust_MX,GAdjust_MM,GAdjust_XX)
jid = lower(JobID);

if contains(jid,{'favrock2'})
    GAdjust_MX = [0 0 1];
    GAdjust_MM = [-3 0.6373 0.005];
    GAdjust_XX = [0 0 1];
elseif contains(jid,{'favrock' 'default' 'def'}) % default
    GAdjust_MX = [0 0 1];
    GAdjust_MM = [0 0 1];
    GAdjust_XX = [0 0 1];
elseif contains(jid,{'favwurtz'})
    GAdjust_MX = [0 0 1];
    GAdjust_MM = [0 0 1];
    GAdjust_XX = [-5 0.7656 0.018];
elseif contains(jid,{'favfiv'})
    GAdjust_MX = [-5 0.5009 0.018];
    GAdjust_MM = [0 0 1];
    GAdjust_XX = [0 0 1];
elseif contains(jid,{'favnias'})
    GAdjust_MX = [0 0 1];
    GAdjust_MM = [-5 0.4769 0.018];
    GAdjust_XX = [0 0 1];
elseif contains(jid,{'favbbeo'})
    GAdjust_MX = [-2 0.5988 0.005];
    GAdjust_MM = [-2 0.4716 0.005];
    GAdjust_XX = [-2 0.4616 0.005];
elseif contains(jid,{'favcscl'})
    GAdjust_MX = [-3 0.2876 0.005; -3 0.5507 0.005];
    GAdjust_MM = [-3 0.3321 0.005];
    GAdjust_XX = [-3 0.3321 0.005];
elseif contains(jid,{'favsphal'})
    GAdjust_MX = [-5 0.4695 0.018];
    GAdjust_MM = [0 0 1];
    GAdjust_XX = [0 0 1];
else
    return
end

end