function LatticeLength = MinLatticeLength(a_vec,b_vec,c_vec)

lattice_vecs = zeros(1,27);
jdx=1;
for ii = [0,-1,1]
    for jj = [0,-1,1]
        for kk = [0,-1,1]
            lattice_vecs(jdx) = norm(ii*a_vec+jj*b_vec+kk*c_vec); % note first entry is always 0
            jdx=jdx+1;
        end
    end
end
LatticeLength = min(lattice_vecs(2:end)); % nm

end