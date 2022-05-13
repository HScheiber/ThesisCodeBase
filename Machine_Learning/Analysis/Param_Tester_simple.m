% Param_Tester_simple
clear;

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};

for idx = 1:length(Salts)
    Salt = Salts{idx};
    Theory = 'TF';
    Fix_Charge = true;
    Additivity = true;
    Additional_MM_Disp = false;
    Additional_GAdjust = {};
    SigmaEpsilon = false;
    C6Damp = Init_C6Damping_Object;
    Additional_Function = Init_Additional_Function;
    Fix_C8 = false;
    Fix_Alpha = false;
    Loss_Options = init_loss_options;

    Param = ones(1,12);
    Model_name = [Salt '_' Theory];
    fname = [Model_name '_reference.mat'];
    resave = true;

    if ispc
        Models_Directory = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
    else
        Models_Directory = '/home/user/project/BO_Models';
    end

    %% Calculation options
    % Note: make sure rocksalt is always included!
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
    Parallel_LiX_Minimizer = true;
    Parallel_Struct_Min    = false;

    tc = tic;
    [Loss,~,Minimization_Data] = LiX_Minimizer(Salt,Structures,...
        Theory,Parallel_LiX_Minimizer,Parallel_Struct_Min,Loss_Options,...
        Fix_Charge,Additivity,Param,'Additional_MM_Disp',Additional_MM_Disp,...
        'GAdjust',Additional_GAdjust,'SigmaEpsilon',SigmaEpsilon,'Verbose',true,...
        'Extra_Properties',true,'C6Damp',C6Damp,'Additional_Function',Additional_Function,...
        'Fix_Alpha',Fix_Alpha,'Fix_C8',Fix_C8);
    t = toc(tc);

    if SigmaEpsilon && Additivity && strcmp(Theory,'JC')        
        % Sigma scaling
        sigma_MM = Param(1);
        sigma_XX = Param(2);

        % Epsilon scaling
        Epsilon_MM = Param(3);
        Epsilon_XX = Param(4);

        Sigma_MX   = (sigma_MM + sigma_XX)/2;
        Epsilon_MX = sqrt(Epsilon_MM*Epsilon_XX);

        Pars(1:2) = Param(1:2);
        Pars(3)   = Sigma_MX;
        Pars(4:5) = Param(3:4);
        Pars(6)   = Epsilon_MX;

        if ~Fix_Charge
            Pars(7) = Param(5);
            if Additional_MM_Disp
                Pars(4) = Pars(4) + Param(6);
            end
        else
            Pars(7) = 1;
            if Additional_MM_Disp
                Pars(4) = Pars(4) + Param(5);
            end
        end


    elseif Additivity && strcmp(Theory,'JC')

        Scale.D.MM = Param(1);
        Scale.D.XX = Param(2);

        % Repulsion
        Scale.R.MM = Param(3);
        Scale.R.XX = Param(4);

        PotSettings = Initialize_MD_Settings;
        PotSettings.Salt = Salt;
        PotSettings.S = Init_Scaling_Object;
        [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

        % Unscaled
        MX_Epsilon = MXParams.epsilon;
        MX_Sigma   = MXParams.sigma;

        MX_R = 4*MX_Epsilon*MX_Sigma^12;
        MX_D = 4*MX_Epsilon*MX_Sigma^6;

        % Scaled
        MM_Epsilon = MMParams.epsilon*(Scale.D.MM^2)*(1/Scale.R.MM);
        MM_Sigma = MMParams.sigma*(1/(Scale.D.MM^(1/6)))*(Scale.R.MM^(1/6));

        XX_Epsilon = XXParams.epsilon*(Scale.D.XX^2)*(1/Scale.R.XX);
        XX_Sigma = XXParams.sigma*(1/(Scale.D.XX^(1/6)))*(Scale.R.XX^(1/6));

        MX_Epsilon = sqrt(MM_Epsilon*XX_Epsilon);
        MX_Sigma   = (MM_Sigma + XX_Sigma)/2;

        MX_R_scaled = 4*MX_Epsilon*MX_Sigma^12;
        MX_D_scaled = 4*MX_Epsilon*MX_Sigma^6;

        Scale.D.MX = MX_D_scaled/MX_D;
        Scale.R.MX = MX_R_scaled/MX_R;

        Pars(1:2) = Param(1:2);
        Pars(3)   = Scale.D.MX;
        Pars(4:5) = Param(3:4);
        Pars(6)   = Scale.R.MX;

        if ~Fix_Charge
            Pars(7) = Param(5);
            if Additional_MM_Disp
                Pars(1) = Pars(1) + Param(6);
            end
        else
            Pars(7) = 1;
            if Additional_MM_Disp
                Pars(1) = Pars(1) + Param(5);
            end
        end

    elseif strcmp(Theory,'JC')
        Pars(1:6) = Param(1:6);
        if Fix_Charge
            Pars(7) = 1;
        else
            Pars(7) = Param(7);
        end
    else % TF model

        Pars(1:3) = Param(1:3);

        if Fix_C8
            Scale.D6D.MM = Pars(1);
            Scale.D6D.XX = Pars(2);
            Scale.D6D.MX = Pars(3);

            % Calculate value of C8 using recursive relations
            [Metal,Halide] = Separate_Metal_Halide(Salt);

            % Default TF params: C in units of kJ/mol nm^6, D in units of kJ/mol nm^8
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Salt;
            PotSettings.S = Init_Scaling_Object;
            [TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(PotSettings);

            % Conversion factors
            Bohr_nm = 0.0529177; % a_0 - > Angstrom
            c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
            J_kJ = 1e-3; % J - > kJ
            Ha_kJmol = 2625.4999; % Ha - > kJ/mol
            c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
            c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

            % Factor used to calculate C8
            sqrt_Q.Li = 5.019869340000000;
            sqrt_Q.Na = 6.585855360000000;
            sqrt_Q.K  = 7.977627530000000;
            sqrt_Q.Rb = 9.554616980000000;
            sqrt_Q.Cs = 11.02204549000000;
            sqrt_Q.F  = 2.388252500000000;
            sqrt_Q.Cl = 3.729323560000000;
            sqrt_Q.Br = 4.590896470000000;
            sqrt_Q.I  = 5.533218150000000;

            % Calculate Scaled C8 using recursion relation from D3 paper
            C8.MM = 3.0*(Scale.D6D.MM*TF_MM.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Metal)*c8units; % in kJ/mol nm^8
            C8.XX = 3.0*(Scale.D6D.XX*TF_XX.C/c6units)*sqrt_Q.(Halide)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8
            C8.MX = 3.0*(Scale.D6D.MX*TF_MX.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8

            % Update the scaling
            Pars(4) = C8.MM/TF_MM.D;
            Pars(5) = C8.XX/TF_XX.D;
            Pars(6) = C8.MX/TF_MX.D;
            pidx = 3;
        else
            Pars(4:6) = Param(4:6);
            pidx = 6;
        end

        % Alpha (TF exponential steepness repulsive parameter)
        if Fix_Alpha
            Pars(7:9) = 1;
        else
            Pars(7) = Param(pidx+1);
            Pars(8) = Param(pidx+2);
            Pars(9) = Param(pidx+3);
            pidx = pidx+3;
        end

        % Repulsive wall prefactor
        Pars(10) = Param(pidx+1);
        Pars(11) = Param(pidx+2);
        Pars(12) = Param(pidx+3);
        pidx = pidx+3;

        if Fix_Charge
            Pars(13) = 1;
        else
            Pars(13) = Param(pidx+1);
        end
    end

    % Add the Gaussians Params
    if ~isempty(Additional_GAdjust)
        N_Pars = length(Pars);
        N_param = length(Param);

        N_adj = length(Additional_GAdjust);
        added_params = N_adj*3;

        par_idx = N_Pars+1;
        param_idx = N_param-added_params+1;

        Pars(par_idx:par_idx+added_params-1) = Param(param_idx:N_param);
    end


    if Additional_Function.MM.N >= 0
        Ex_fun = 1;
    else
        Ex_fun = 0;
    end

    if Additional_Function.XX.N >= 0
        Ex_fun = Ex_fun + 1;
    end
    if Additional_Function.MX.N >= 0
        Ex_fun = Ex_fun + 1;
    end

    if Ex_fun > 0
        L_pars = length(Pars);
        Pars(L_pars+1:L_pars+Ex_fun) = Param(end+1-Ex_fun:end);
    end


    En = zeros(size(Structures)); % 'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
    for idx = 1:length(Structures)
        En(idx) = Minimization_Data{idx}.E;
        switch Structures{idx}
            case {'Rocksalt' 'Sphalerite' 'CsCl'}
                En(end+1) = Minimization_Data{idx}.a;
            otherwise
                En(end+1) = Minimization_Data{idx}.a;
                En(end+1) = Minimization_Data{idx}.c;
        end
        disp(Structures{idx})
        disp(Minimization_Data{idx})
    end

    disp(['Current model: ' Model_name])
    disp(['Total Loss: ' num2str(Loss,'%.10f')])
    disp(['Total Time: ' num2str(t,'%.2f') ' sec.'])
    if SigmaEpsilon
        disp(['Parameters (Sigma/Epsilon): ' num2str(Pars,'\t%.15f')])
    else
        disp(['Parameters (Dispersion/Repulsion Scale): ' num2str(Pars,'\t%.15f')])
    end
    disp(['Energies (kJ/mol): ' num2str(En,'\t%.15f')])

    if resave
        vars.Minimization_Data = Minimization_Data;
        vars.loss = Loss;

        saveloc = fullfile(Models_Directory,fname);
        save(saveloc,'-struct','vars')
    end
end