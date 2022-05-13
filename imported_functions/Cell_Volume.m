function Volume = Cell_Volume(a_vec,b_vec,c_vec)

Volume = abs(dot(cross(a_vec,b_vec),c_vec));


end