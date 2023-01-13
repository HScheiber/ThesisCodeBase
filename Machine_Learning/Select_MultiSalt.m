function [Salts,MultiMetal,MultiHalide] = Select_MultiSalt(Settings)
    
    [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
    
    if strcmp(Metal,'M')
        Metals = {'Li' 'Na' 'K' 'Rb' 'Cs'};
        MultiMetal = true;
    else
        Metals = {Metal};
        MultiMetal = false;
    end
    if strcmp(Halide,'X')
        Halides = {'F' 'Cl' 'Br' 'I'};
        MultiHalide = true;
    else
        Halides = {Halide};
        MultiHalide = false;
    end
    
    A = 1:numel(Metals);
    B = 1:numel(Halides);
    [m,n] = ndgrid(A,B);
    Z = [m(:),n(:)];
    S = horzcat(Metals(:,Z(:,1))',Halides(:,Z(:,2))');
    Salts = join(S,'',2);
end