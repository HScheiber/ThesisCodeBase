function p_name = par_name_map(p_name)
    p_name = strrep(p_name,'epsilon_','$\epsilon_{');
    p_name = strrep(p_name,'sigma_','$\sigma_{');
    p_name = strrep(p_name,'gamma_','$\gamma_{');
    p_name = strrep(p_name,'r0_','$R_{0,');
    p_name = strrep(p_name,'M','\textrm{Li}^{+}');
    p_name = strrep(p_name,'X','\textrm{X}^{-}');
    p_name = strrep(p_name,'F','\textrm{F}^{-}');
    p_name = strrep(p_name,'Cl','\textrm{Cl}^{-}');
    p_name = strrep(p_name,'Br','\textrm{Br}^{-}');
    p_name = strrep(p_name,'I','\textrm{I}^{-}');
    p_name = [p_name '}$'];    
end