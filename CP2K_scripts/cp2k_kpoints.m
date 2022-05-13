function [KPOINT_OUT,nrep] = cp2k_kpoints(kpoints,HFX,Structure,ABC,cutoff_radius)

    if HFX
        KPOINT_OUT = ''; % Calculations for HFX only at gamma point
        
        ABC_Split = split(ABC,'   ');
        a = str2double(ABC_Split{1});
        c = str2double(ABC_Split{3});
        
        switch lower(Structure)
            case 'rocksalt'
                min_dim = a*sind(60)/2;
                
                n = ceil(cutoff_radius/min_dim);
                nrep.a = n;
                nrep.b = n;
                nrep.c = n;
                
            case 'wurtzite'
                min_dim_ab = a*sind(60)/2;
                min_dim_c = c/2;
                
                n_ab = ceil(cutoff_radius/min_dim_ab);
                n_c = ceil(cutoff_radius/min_dim_c);
                nrep.a = n_ab;
                nrep.b = n_ab;
                nrep.c = n_c;
                
            case 'fivefive'
                min_dim_ab = a*sind(60)/2;
                min_dim_c = c/2;
                
                n_ab = ceil(cutoff_radius/min_dim_ab);
                n_c = ceil(cutoff_radius/min_dim_c);
                nrep.a = n_ab;
                nrep.b = n_ab;
                nrep.c = n_c;
                
            case 'cscl'
                min_dim = a*sind(60)/2;
                
                n = ceil(cutoff_radius/min_dim);
                nrep.a = n;
                nrep.b = n;
                nrep.c = n;
                
            case 'betabeo'
                min_dim_ab = a/2;
                min_dim_c = c/2;
                
                n_ab = ceil(cutoff_radius/min_dim_ab);
                n_c = ceil(cutoff_radius/min_dim_c);
                nrep.a = n_ab;
                nrep.b = n_ab;
                nrep.c = n_c;
                
            case 'sphalerite'
                min_dim = a*sind(60)/2;
                
                n = ceil(cutoff_radius/min_dim);
                nrep.a = n;
                nrep.b = n;
                nrep.c = n;
                
            case {'nias' 'antinias'}
                min_dim_ab = a*sind(60)/2;
                min_dim_c = c/2;
                
                n_ab = ceil(cutoff_radius/min_dim_ab);
                n_c = ceil(cutoff_radius/min_dim_c);
                nrep.a = n_ab;
                nrep.b = n_ab;
                nrep.c = n_c;
                
        end
        
    else
        kpoint_txt = num2str(kpoints);
        KPOINT_OUT =   ['	&KPOINTS' newline ...
                        '		SCHEME MONKHORST-PACK ' kpoint_txt ' ' kpoint_txt ' ' kpoint_txt newline ...
                        '	&END KPOINTS'];
        nrep.a = 1;
        nrep.b = 1;
        nrep.c = 1;
    end
    
nrep.a = num2str(nrep.a);
nrep.b = num2str(nrep.b);
nrep.c = num2str(nrep.c);

end