function [Li,Na,K,Rb,F,Cl,Br,I] = Prune_Scale_pob_TZVP(home,Prune_M1,Prune_M2,...
    Prune_H1,Prune_H2,Prune_H3,Scale_M1,Scale_M2,Scale_H1,Scale_H2,Scale_H3,...
    Add_M,Add_H,NumM,NumH,current_crystal_type,BSSE_id,Initialize_As_Ions)

load([home filesep 'basis_sets' filesep 'pob_TZVP_CRYSTAL17.mat']) %#ok<LOAD>
RegS = ' ([1-9]|10)\.0 '; % Regex string to replace charges

% Modify Charges if needed
if strcmp(current_crystal_type,'Ion')
    Metal_Charge = Metal_Charge-1; %#ok<*NODEF>
    Halide_Charge = Halide_Charge+1;
elseif strcmp(current_crystal_type,'Atom')
    % Do not modify charges
elseif Initialize_As_Ions
    Metal_Charge = Metal_Charge-1;
    Halide_Charge = Halide_Charge+1;
end

% Replace Charges back into text
Li = strrep(Li,'##CHG##',num2str(Metal_Charge,'%4.1f'));
Na = strrep(Na,'##CHG##',num2str(Metal_Charge,'%4.1f'));
K = strrep(K,'##CHG##',num2str(Metal_Charge,'%4.1f'));
Rb = strrep(Rb,'##CHG##',num2str(Metal_Charge,'%4.1f'));
F = strrep(F,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Cl = strrep(Cl,'##CHG##',num2str(Halide_Charge,'%4.1f'));
Br = strrep(Br,'##CHG##',num2str(Halide_Charge,'%4.1f'));
I = strrep(I,'##CHG##',num2str(Halide_Charge,'%4.1f'));

% Add Extra basis function(s)
Li = strrep(Li,[newline '##EXT##'],Add_M);
Li_NUM = Li_NUM + NumM;

Na = strrep(Na,[newline '##EXT##'],Add_M);
Na_NUM = Na_NUM + NumM;

K = strrep(K,[newline '##EXT##'],Add_M);
K_NUM = K_NUM + NumM;

Rb = strrep(Rb,[newline '##EXT##'],Add_M);
Rb_NUM = Rb_NUM + NumM;

F = strrep(F,[newline '##EXT##'],Add_H);
F_NUM = F_NUM + NumH;

Cl = strrep(Cl,[newline '##EXT##'],Add_H);
Cl_NUM = Cl_NUM + NumH;

Br = strrep(Br,[newline '##EXT##'],Add_H);
Br_NUM = Br_NUM + NumH;

I = strrep(I,[newline '##EXT##'],Add_H);
I_NUM = I_NUM + NumH;

% Prune Metal basis if requested
if Prune_M1
    Li = strrep(Li,[newline '##P1##'],'');
    Li_NUM = Li_NUM - 1;
    
    Na = strrep(Na,[newline '##P1##'],'');
    Na_NUM = Na_NUM - 1;
    
    K = strrep(K,[newline '##P1##'],'');
    K_NUM = K_NUM - 1;
    
    Rb = strrep(Rb,[newline '##P1##'],'');
    Rb_NUM = Rb_NUM - 1;
    
else
    Li = strrep(Li,'##P1##',Li_P1);
    
    Na = strrep(Na,'##P1##',Na_P1);
    
    K = strrep(K,'##P1##',K_P1);
    
    Rb = strrep(Rb,'##P1##',Rb_P1);
end

if Prune_M2
    Li = strrep(Li,[newline '##P2##'],'');
    Li_NUM = Li_NUM - 1;
    
    Na = strrep(Na,[newline '##P2##'],'');
    Na_NUM = Na_NUM - 1;
    
    K = strrep(K,[newline '##P2##'],'');
    K_NUM = K_NUM - 1;
    
    Rb = strrep(Rb,[newline '##P2##'],'');
    Rb_NUM = Rb_NUM - 1;
else
    Li = strrep(Li,'##P2##',Li_P2);
    
    Na = strrep(Na,'##P2##',Na_P2);
    
    K = strrep(K,'##P2##',K_P2);
    
    Rb = strrep(Rb,'##P2##',Rb_P2);
end

Li = strrep(Li,'##NUM##',num2str(Li_NUM));
Na = strrep(Na,'##NUM##',num2str(Na_NUM));
K = strrep(K,'##NUM##',num2str(K_NUM));
Rb = strrep(Rb,'##NUM##',num2str(Rb_NUM));

% Prune Halide basis if requested
if Prune_H1
    F = strrep(F,[newline '##P1##'],'');
    Cl = strrep(Cl,[newline '##P1##'],'');
    Br = strrep(Br,[newline '##P1##'],'');
    I = strrep(I,[newline '##P1##'],'');
    
    F_NUM = F_NUM - 1;
    Cl_NUM = Cl_NUM - 1;
    Br_NUM = Br_NUM - 1;
    I_NUM = I_NUM - 1;
else
    F = strrep(F,'##P1##',F_P1);
    Cl = strrep(Cl,'##P1##',Cl_P1);
    Br = strrep(Br,'##P1##',Br_P1);
    I = strrep(I,'##P1##',I_P1);
end

if Prune_H2
    F = strrep(F,[newline '##P2##'],'');
    Cl = strrep(Cl,[newline '##P2##'],'');
    Br = strrep(Br,[newline '##P2##'],'');
    I = strrep(I,[newline '##P2##'],'');
    
    F_NUM = F_NUM - 1;
    Cl_NUM = Cl_NUM - 1;
    Br_NUM = Br_NUM - 1;
    I_NUM = I_NUM - 1;
else
    F = strrep(F,'##P2##',F_P2);
    Cl = strrep(Cl,'##P2##',Cl_P2);
    Br = strrep(Br,'##P2##',Br_P2);
    I = strrep(I,'##P2##',I_P2);
end

if Prune_H3
    F = strrep(F,[newline '##P3##'],'');
    Cl = strrep(Cl,[newline '##P3##'],'');
    Br = strrep(Br,[newline '##P3##'],'');
    I = strrep(I,[newline '##P3##'],'');
    
    F_NUM = F_NUM - 1;
    Cl_NUM = Cl_NUM - 1;
    Br_NUM = Br_NUM - 1;
    I_NUM = I_NUM - 1;
else
    F = strrep(F,'##P3##',F_P3);
    Cl = strrep(Cl,'##P3##',Cl_P3);
    Br = strrep(Br,'##P3##',Br_P3);
    I = strrep(I,'##P3##',I_P3);
end
F = strrep(F,'##NUM##',num2str(F_NUM));
Cl = strrep(Cl,'##NUM##',num2str(Cl_NUM));
Br = strrep(Br,'##NUM##',num2str(Br_NUM));
I = strrep(I,'##NUM##',num2str(I_NUM));

% Scale basis set
Li_exponent1 = num2str(Li_X1*Scale_M1,'%1.10f');
Li_exponent2 = num2str(Li_X2*Scale_M2,'%1.10f');
Li = strrep(Li,'##X1##',Li_exponent1);
Li = strrep(Li,'##X2##',Li_exponent2);

Na_exponent1 = num2str(Na_X1*Scale_M1,'%1.10f');
Na_exponent2 = num2str(Na_X2*Scale_M2,'%1.10f');
Na = strrep(Na,'##X1##',Na_exponent1);
Na = strrep(Na,'##X2##',Na_exponent2);

K_exponent1 = num2str(K_X1*Scale_M1,'%1.10f');
K_exponent2 = num2str(K_X2*Scale_M2,'%1.10f');
K = strrep(K,'##X1##',K_exponent1);
K = strrep(K,'##X2##',K_exponent2);

Rb_exponent1 = num2str(Rb_X1*Scale_M1,'%1.10f');
Rb_exponent2 = num2str(Rb_X2*Scale_M2,'%1.10f');
Rb = strrep(Rb,'##X1##',Rb_exponent1);
Rb = strrep(Rb,'##X2##',Rb_exponent2);

F_exponent1 = num2str(F_X1*Scale_H1,'%1.10f');
F_exponent2 = num2str(F_X2*Scale_H2,'%1.10f');
F_exponent3 = num2str(F_X3*Scale_H3,'%1.10f');
F = strrep(F,'##X1##',F_exponent1);
F = strrep(F,'##X2##',F_exponent2);
F = strrep(F,'##X3##',F_exponent3);

Cl_exponent1 = num2str(Cl_X1*Scale_H1,'%1.10f');
Cl_exponent2 = num2str(Cl_X2*Scale_H2,'%1.10f');
Cl_exponent3 = num2str(Cl_X3*Scale_H3,'%1.10f');
Cl = strrep(Cl,'##X1##',Cl_exponent1);
Cl = strrep(Cl,'##X2##',Cl_exponent2);
Cl = strrep(Cl,'##X3##',Cl_exponent3);

Br_exponent1 = num2str(Br_X1*Scale_H1,'%1.10f');
Br_exponent2 = num2str(Br_X2*Scale_H2,'%1.10f');
Br_exponent3 = num2str(Br_X3*Scale_H3,'%1.10f');
Br = strrep(Br,'##X1##',Br_exponent1);
Br = strrep(Br,'##X2##',Br_exponent2);
Br = strrep(Br,'##X3##',Br_exponent3);

I_exponent1 = num2str(I_X1*Scale_H1,'%1.10f');
I_exponent2 = num2str(I_X2*Scale_H2,'%1.10f');
I_exponent3 = num2str(I_X3*Scale_H3,'%1.10f');
I = strrep(I,'##X1##',I_exponent1);
I = strrep(I,'##X2##',I_exponent2);
I = strrep(I,'##X3##',I_exponent3);

% For custom pair BSSE
if strcmp(BSSE_id,'Halide') && strcmp(current_crystal_type,'Pair')
    Li_ATNUM = 0;
    Na_ATNUM = 0;
    K_ATNUM = 0;
    Rb_ATNUM = 0;
    Rb = strrep(Rb,['##INPUT##' newline],'');
    
    % Remove electrons
    Li = regexprep(Li,RegS,' 0.0 ');
    Na = regexprep(Na,RegS,' 0.0 ');
    K = regexprep(K,RegS,' 0.0 ');
    Rb = regexprep(Rb,'(0 [0-3] [1-9] )([1-9]|10)(\.*0*)( 1.0)','$10$4');
elseif strcmp(BSSE_id,'Metal') && strcmp(current_crystal_type,'Pair')
    F_ATNUM = 0;
    Cl_ATNUM = 0;
    Br_ATNUM = 0;
    I_ATNUM = 0;
    I = strrep(I,['##INPUT##' newline],'');
    
    % Remove electrons
    F = regexprep(F,RegS,' 0.0 ');
    Cl = regexprep(Cl,RegS,' 0.0 ');
    Br = regexprep(Br,RegS,' 0.0 ');
    I = regexprep(I,'(0 [0-3] [1-9] )([1-9]|10)(\.*0*)( 1.0)','$10$4');
end

Li = strrep(Li,'##ATNUM##',num2str(Li_ATNUM,'%u'));
Na = strrep(Na,'##ATNUM##',num2str(Na_ATNUM,'%u'));
K = strrep(K,'##ATNUM##',num2str(K_ATNUM,'%u'));
Rb = strrep(Rb,'##ATNUM##',num2str(Rb_ATNUM,'%u'));
Rb = strrep(Rb,'##INPUT##',Rb_INPUT);
F = strrep(F,'##ATNUM##',num2str(F_ATNUM,'%u'));
Cl = strrep(Cl,'##ATNUM##',num2str(Cl_ATNUM,'%u'));
Br = strrep(Br,'##ATNUM##',num2str(Br_ATNUM,'%u'));
I = strrep(I,'##ATNUM##',num2str(I_ATNUM,'%u'));
I = strrep(I,'##INPUT##',I_INPUT);

end