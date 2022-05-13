salt_types = "Li" + ["Cl" "Br" "I"];
crystal_types = ["Wurtzite" "Rocksalt"];

for salt=salt_types
    mkdir(salt)
    cd(salt)
    for crys=crystal_types
        mkdir(crys)
        cd(crys)
        F_submit_crystal_PES_Array(salt,crys)
        cd("..")
    end
    cd("..")
end
