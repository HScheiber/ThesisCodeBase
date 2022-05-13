function [RDF,CRDF] = calc_RDF(wb,wsl,ech,gmx,Traj_file,Trajectory_Out_file,Tpr_file,Start_Point,End_Point,...
    TimePerFrame,pbc_mod,pbc_on,RDF_Out_file,ref,sel,CRDF_Out_file)

    if pbc_on
        Traj_in_file = Traj_file;
    else
        
        % Grab trajectory data from results
        gmx_command = [wsl ech gmx ' trjconv -f ' windows2unix(Traj_file)...
                    ' -o ' windows2unix(Trajectory_Out_file) ' -s ' windows2unix(Tpr_file) ...
                    ' -b ' Start_Point ' -e ' End_Point ' -dt ' TimePerFrame pbc_mod];
                
        % If PBC is off, need to convert trajectory first
        err = system(gmx_command);
        if err ~= 0
            warndlg('Failed to collect trajectory for RDF analysis with PBC off.')
            return
        end
        Traj_in_file = Trajectory_Out_file;
    end
    
    % Check if cluster
    [~,traj_name,~] = fileparts(Traj_file);
    if ~isempty(regexp(traj_name,'_N[0-9]+','ONCE'))
        norm = 'none';
    else
        norm = 'rdf';
    end

    % Grab rdf data from results
    gmx_command = [wsl gmx ' rdf -f ' windows2unix(Traj_in_file)...
        ' -o ' windows2unix(RDF_Out_file) ' -s ' windows2unix(Tpr_file) ...
        ' -b ' Start_Point ' -e ' End_Point ' -rmax 2' ...
        ' -ref ' ref ' -sel ' sel ' -dt ' TimePerFrame ' -cn ' windows2unix(CRDF_Out_file) ...
        ' -norm ' norm];

    err = system(gmx_command);
    if err ~= 0
        RDF = nan(1,2);
        CRDF = nan(1,2);
        return
    end

    RDF = import_xvg(RDF_Out_file); % Gather RDF
    CRDF = import_xvg(CRDF_Out_file); % Gather CRDF

    % Remove the temporary output files
    delete(RDF_Out_file)
    delete(CRDF_Out_file)
    increment(wb)
end