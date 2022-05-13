function Energy = GrabEnergy(Directory,FileBase,varargin)

    p = inputParser;
    p.FunctionName = 'GrabEnergy';
    addOptional(p,'Extra_Properties',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
    parse(p,varargin{:});
    Extra_Properties = p.Results.Extra_Properties;

    % Find xvg file
    enfile_info = dir(fullfile(Directory,...
        [FileBase 'energies.xvg']));

    % Find total number of atoms .mat file
    totatoms_info = dir(fullfile(Directory,...
        [FileBase '.mat']));

    if isempty(enfile_info) || isempty(totatoms_info)
        error('Optimization step failed to produce required output files.');
    end

    % Get total atoms filename
    Atoms_file = fullfile(Directory,totatoms_info.name);

    % Open the mat file to get total number of atoms/molecules
    N_info = load(Atoms_file,'-mat','N_Cell','N_total');
    N_total_mols = N_info.N_total/2;

    % Get energy filename
    Energy_file = fullfile(Directory,enfile_info.name);

    % Import energy file as text
    Energy_text = fileread(Energy_file);

    % Find the first energy line
    Energies_Initial = regexp(Energy_text,'    0.000000.+?\n','match');
    
    if Extra_Properties
        Energy_list = textscan(Energies_Initial{1},'%*f %f %f %f %f %f %f %f %f %f %f',...
            'Delimiter',' ','MultipleDelimsAsOne',true);
        Energy_array = cell2mat(Energy_list)./N_total_mols;

        Energy.LJSR = Energy_array(1);
        Energy.CoulSR = Energy_array(2);
        Energy.CoulLR = Energy_array(3);
        Energy.Pot = Energy_array(4);
        Energy.CoulSR_MM = Energy_array(5);
        Energy.LJSR_MM = Energy_array(6);
        Energy.CoulSR_MX = Energy_array(7);
        Energy.LJSR_MX = Energy_array(8);
        Energy.CoulSR_XX = Energy_array(9);
        Energy.LJSR_XX = Energy_array(10);
    else
        Energy_list = textscan(Energies_Initial{1},'%*f %f',...
            'Delimiter',' ','MultipleDelimsAsOne',true);
        Energy = cell2mat(Energy_list)./N_total_mols;
    end
end