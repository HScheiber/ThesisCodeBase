Labels = ["Cl" "Br" "I"];
for Label=Labels
    mkdir(Label + "_Ion")
    cd(Label + "_Ion")
    submit_benchmarking_theory_basis(Label)
    cd("..")
end
