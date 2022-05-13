function I = Prune_Scale_ECP28MDF_VQZ(home,Add_H,NumH,...
    current_crystal_type,BSSE_id,Initialize_As_Ions)

load([home filesep 'basis_sets' filesep 'ECP28MDF_VQZ_CRYSTAL17.mat']) %#ok<LOAD>
RegS = ' ([1-9]|10)\.0 '; % Regex string to replace charges

% Modify Charges if needed
if strcmp(current_crystal_type,'Ion')
    Halide_Charge = Halide_Charge+1;%#ok<*NODEF>
elseif strcmp(current_crystal_type,'Atom')
    % Do not modify charges
elseif Initialize_As_Ions
    Halide_Charge = Halide_Charge+1;
end

% Replace Charges back into text
I = strrep(I,'##CHG##',num2str(Halide_Charge,'%4.1f'));

% Add Extra basis function(s)
I = strrep(I,[newline '##EXT##'],Add_H);
I_NUM = I_NUM + NumH;

I = strrep(I,'##NUM##',num2str(I_NUM));

% Scale basis set

% For custom pair BSSE
if strcmp(BSSE_id,'Halide') && strcmp(current_crystal_type,'Pair')

elseif strcmp(BSSE_id,'Metal') && strcmp(current_crystal_type,'Pair')
    I_ATNUM = 0;
    I = strrep(I,['##INPUT##' newline],'');
    
    % Remove electrons
    I = regexprep(I,'(0 [0-3] [1-9] )([1-9]|10)(\.*0*)( 1.0)','$10$4');
end
I = strrep(I,'##ATNUM##',num2str(I_ATNUM,'%u'));
I = strrep(I,'##INPUT##',I_INPUT);

end