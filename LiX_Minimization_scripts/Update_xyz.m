function Geometry = Update_xyz(Geometry)

    a_vec = [Geometry.a/10 0 0]*Geometry.Transform;
    b_vec = [0 Geometry.b/10 0]*Geometry.Transform;
    c_vec = [0 0 Geometry.c/10]*Geometry.Transform;
    TM = [a_vec ; b_vec ; c_vec];
    
    % Update geometry
    Geometry.xyz = Geometry.FC*TM;
    
    % Update boxcoords
    Geometry.boxcoords([1,4,5]) = num2cell(a_vec);
    Geometry.boxcoords([6,2,7]) = num2cell(b_vec);
    Geometry.boxcoords([8,9,3]) = num2cell(c_vec);
    
end