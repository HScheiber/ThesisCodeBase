function Minimization_Data = Initialize_Minimization_Data(Structures,Salt,Theory)

N = length(Structures);
Minimization_Data = cell(N,1);
for idx = 1:N
    this.LJSR = nan;
    this.CoulSR = nan;
    this.CoulLR = nan;
    this.Pot = nan;
    this.CoulSR_MM = nan;
    this.LJSR_MM = nan;
    this.CoulSR_MX = nan;
    this.LJSR_MX = nan;
    this.CoulSR_XX = nan;
    this.LJSR_XX = nan;
    this.N_atoms = nan;
    this.a = nan;
    this.b = nan;
    this.c = nan;
    this.FC_Metal = [nan nan nan];
    this.FC_Halide = [nan nan nan];
    this.V = nan;
    this.E = nan;
    this.Salt = Salt;
    this.Structure = Structures{idx};
    this.Model = Theory;
    Minimization_Data{idx} = this;
end

end