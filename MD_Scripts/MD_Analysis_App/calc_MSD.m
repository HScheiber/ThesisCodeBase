function [MSD,DiffCoeff,DiffCoeff_e] = calc_MSD(wb,wsl,ech,gmx,Traj_file,Trajectory_Out_file,Tpr_file,Start_Point,End_Point,...
    TimePerFrame,pbc_mod,pbc_on,MSD_Out_file,ref)
   
    if pbc_on
        Traj_in_file = Traj_file;
    else
        % Grab trajectory data from results
        gmx_command = [wsl ech gmx ' trjconv -f ' windows2unix(Traj_file)...
                    ' -o ' windows2unix(Trajectory_Out_file) ' -s ' windows2unix(Tpr_file) ...
                    ' -b ' Start_Point ' -e ' End_Point ' -dt ' TimePerFrame pbc_mod];
                
        % If PBC is off, need to convert trajectory first?
        err = system(gmx_command);
        if err ~= 0
            warndlg('Failed to collect trajectory for MSD analysis with PBC off.')
            return
        end
        Traj_in_file = Trajectory_Out_file;
    end

    % Grab msd data from results
    gmx_command = [wsl ref gmx ' msd -f ' windows2unix(Traj_in_file)...
        ' -o ' windows2unix(MSD_Out_file) ' -s ' windows2unix(Tpr_file) ...
        ' -b ' Start_Point ' -e ' End_Point ...
        ' -trestart ' TimePerFrame ];

    [err,output] = system(gmx_command);
    if err ~= 0
        MSD = nan(1,2);
        return
    end
    
    % Grab diffusion coefficient and error
    DiffCoeffs = regexp(output,'D\[.+?\] *(([0-9]|\.)+) *\(\+/- *(([0-9]|\.)+)\) *(([0-9]|\.|e|-)+)','tokens','once');
    DiffCoeff_e = roundsd(str2double(DiffCoeffs{2})*str2double(DiffCoeffs{3})*1e5,1); % units nm^2/ns
    N = -(floor(log10(abs(DiffCoeff_e))));
    DiffCoeff = round(str2double(DiffCoeffs{1})*str2double(DiffCoeffs{3})*1e5,N); % units nm^2/ns
    
    MSD = import_xvg(MSD_Out_file); % Gather MSD

    % Remove the temporary output files
    delete(MSD_Out_file)
    increment(wb)
end