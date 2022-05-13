function CRDamping = Init_CRDamping_Object

    CRDamping.MX.r_d = -1; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
    CRDamping.MX.b = -1; % sigmoid "steepness" for damping
    CRDamping.MM.r_d = -1; % LiCl = 0.21 , LiI = 0.24
    CRDamping.MM.b  = -1; % 75
    CRDamping.XX.r_d = -1; 
    CRDamping.XX.b  = -1;
    
end