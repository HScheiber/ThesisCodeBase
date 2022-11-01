function [Topology_Template,MDP_Template] = Polarize_Inputs(Settings,Topology_Template,MDP_Template)
    
    [Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
    
    %% Core and Shell charges
    % Net charges
    q.M =  Settings.S.Q; % atomic
    q.X  = -Settings.S.Q; % atomic

    q.core.M = Settings.S.QcoreM;
    q.shell.M = q.M - Settings.S.QcoreM;

    q.core.X = Settings.S.QcoreX;
    q.shell.X = q.X - Settings.S.QcoreX;
    
    alpha.M = Settings.S.PM;
    alpha.X = Settings.S.PX;
    
    %% Deal with energy groups in MDP file
    energy_groups = [Metal ' ' Metal '_s' ' ' Halide ' ' Halide '_s'];
    MDP_Template = strrep(MDP_Template,'##ENERGYGRPS##',energy_groups);
    MDP_Template = regexprep(MDP_Template,'nstcalcenergy * = *[0-9]+ *','nstcalcenergy            = 1                 ');
    MDP_Template = [MDP_Template newline newline...
        '; maximum number of iterations for optimizing the shell positions and the flexible constraints' newline ...
        'niter                    = ' num2str(Settings.niter_polarization)];
    MDP_Template = [MDP_Template newline newline...
        '; [kJ mol-1 nm-1] the minimization is converged when the maximum force is smaller than this value' newline ...
        'emtol                    = ' num2str(Settings.emtol_polarization,'%.2E')];
    
    %% Deal with topology
    
    % atomtypes
    Topology_Template = regexprep(Topology_Template,'\r','');
    atomtypes = regexp(Topology_Template,'\[ *atomtypes *\].+?\n{2}','match','once');
    metal_info = regexp(atomtypes,['  ' Metal ' *' Metal '.+?\n'],'match','once');
    metal_info = [metal_info '  ' pad([Metal '_s'],6) 'X          0          0.000   ' ...
        num2str(q.shell.M,'%+.3f') '       S      0.0              0.0' newline];
    
    halide_info = regexp(atomtypes,['  ' Halide ' *' Halide '.+?\n'],'match','once');
    halide_info = [halide_info '  ' pad([Halide '_s'],6) 'X          0          0.000   ' ...
        num2str(q.shell.X,'%+.3f') '       S      0.0              0.0' newline];
    
    atomtypes = regexprep(atomtypes,['  ' Metal ' *' Metal '.+?\n'],metal_info);
    atomtypes = regexprep(atomtypes,['  ' Halide ' *' Halide '.+?\n'],halide_info);
    Topology_Template = regexprep(Topology_Template,'\[ *atomtypes *\].+?\n{2}',atomtypes);
    
    % Metal Molecule and polarization
    moleculetypes = regexp(Topology_Template,['\[ *moleculetype *\].+?' Metal ' +' Metal '.+?\n{2}'],'match','once');
    metal_info = regexp(moleculetypes,['1 *' Metal '.+?\n'],'match','once');
    metal_info = [metal_info '2           ' pad([Metal '_s'],11) '1          ' pad(Metal,11) pad([Metal '_s'],11) ...
        '2         ' pad(num2str(q.shell.M,'%+.3f'),11) newline newline ...
        '[ polarization ]' newline ';atom shell type alpha' newline '1     2     1    ' ...
        num2str(alpha.M,'%.4E') newline];
    moleculetypes = regexprep(moleculetypes,['1 *' Metal '.+?\n'],metal_info);
    Topology_Template = regexprep(Topology_Template,['\[ *moleculetype *\].+?' Metal ' +' Metal '.+?\n{2}'],moleculetypes);
    
    % Halide molecule and polarization
    moleculetypes = regexp(Topology_Template,['\[ *moleculetype *\].+?' Halide ' +' Halide '.+?\n{2}'],'match','once');
    halide_info = regexp(moleculetypes,['1 *' Halide '.+?\n'],'match','once');
    halide_info = [halide_info '2           ' pad([Halide '_s'],11) '1          ' pad(Halide,11) pad([Halide '_s'],11) ...
        '2         ' pad(num2str(q.shell.X,'%+.3f'),11) newline newline ...
        '[ polarization ]' newline ';atom shell type alpha' newline '1     2     1    ' ...
        num2str(alpha.X,'%.4E') newline];
    moleculetypes = regexprep(moleculetypes,['1 *' Halide '.+?\n'],halide_info);
    Topology_Template = regexprep(Topology_Template,['\[ *moleculetype *\].+?' Halide ' +' Halide '.+?\n{2}'],moleculetypes);
    
end