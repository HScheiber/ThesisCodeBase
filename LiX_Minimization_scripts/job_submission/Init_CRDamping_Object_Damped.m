function CRDamping = Init_CRDamping_Object_Damped

    CRDamping.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
    CRDamping.MX.b = 100; % sigmoid "steepness" for damping
    CRDamping.MM.r_d = 0.15; % LiCl = 0.21 , LiI = 0.24
    CRDamping.MM.b  = 100; % 75
    CRDamping.XX.r_d = 0.15; 
    CRDamping.XX.b  = 100;
    
end