function params = bayesopt_params(Model)
    

    switch Model.Theory
    case 'TF'
        if ~isfield(Model,'InnerRange')
            Model.InnerRange = true;
        end
        SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
        
        if Model.SigmaEpsilon
            if Model.InnerRange % For TF model, inner range is when gamma is [0, 48/7]
                % r_0
                Sr0MM = optimizableVariable('r0_MM',[0 1],'Type','real'); % Units: nm
                Sr0XX = optimizableVariable('r0_XX',[0 1],'Type','real'); % Units: nm
                Sr0MX = optimizableVariable('r0_MX',[0 1],'Type','real'); % Units: nm
                
                % epsilon
                SepsilonMM = optimizableVariable('epsilon_MM',[1e-4 1e12],'Type','real','Transform','log'); % Units: kJ/mol
                SepsilonXX = optimizableVariable('epsilon_XX',[1e-4 1e12],'Type','real','Transform','log'); % Units: kJ/mol
                SepsilonMX = optimizableVariable('epsilon_MX',[1e-4 1e12],'Type','real','Transform','log'); % Units: kJ/mol
                
                % gamma
                SgammaMM = optimizableVariable('gamma_MM',[0.09 6],'Type','real'); % Units: kJ/mol
                SgammaXX = optimizableVariable('gamma_XX',[0.09 6],'Type','real'); % Units: kJ/mol
                SgammaMX = optimizableVariable('gamma_MX',[0.09 6],'Type','real'); % Units: kJ/mol
            
            else
                % r_0
                Sr0MM = optimizableVariable('r0_MM',[0 1],'Type','real'); % Units: nm
                Sr0XX = optimizableVariable('r0_XX',[0 1],'Type','real'); % Units: nm
                Sr0MX = optimizableVariable('r0_MX',[0 1],'Type','real'); % Units: nm
                
                % epsilon
                SepsilonMM = optimizableVariable('epsilon_MM',[0 400],'Type','real'); % Units: kJ/mol
                SepsilonXX = optimizableVariable('epsilon_XX',[0 400],'Type','real'); % Units: kJ/mol
                SepsilonMX = optimizableVariable('epsilon_MX',[0 400],'Type','real'); % Units: kJ/mol
                
                % gamma
                SgammaMM = optimizableVariable('gamma_MM',[48/7 30],'Type','real'); % Units: kJ/mol
                SgammaXX = optimizableVariable('gamma_XX',[48/7 30],'Type','real'); % Units: kJ/mol
                SgammaMX = optimizableVariable('gamma_MX',[48/7 30],'Type','real'); % Units: kJ/mol
            end
            
            
            if Model.Additivity
                params = [Sr0MM,Sr0XX,SepsilonMM,SepsilonXX,SgammaMX];
            else
                params = [Sr0MM,Sr0XX,Sr0MX,SepsilonMM,SepsilonXX,SepsilonMX,SgammaMM,SgammaXX,SgammaMX];
            end
        else
            % C6 dispersion parameters
            SD6MM = optimizableVariable('SD6MM',[0 1000],'Type','real');
            SD6XX = optimizableVariable('SD6XX',[0 1000],'Type','real');
            SD6MX = optimizableVariable('SD6MX',[0 1000],'Type','real');
            
            % C8 dispersion parameters
            SD8MM = optimizableVariable('SD8MM',[0 1000],'Type','real');
            SD8XX = optimizableVariable('SD8XX',[0 1000],'Type','real');
            SD8MX = optimizableVariable('SD8MX',[0 1000],'Type','real');
            
            % Alpha exponential parameters
            SAMM = optimizableVariable('SAMM',[0.1 10],'Type','real');
            SAXX = optimizableVariable('SAXX',[0.1 10],'Type','real');
            SAMX = optimizableVariable('SAMX',[0.1 10],'Type','real');
            
            % Repulsive wall prefactors
            SRMM = optimizableVariable('SRMM',[0.1 50],'Type','real');
            SRXX = optimizableVariable('SRXX',[0.1 50],'Type','real');
            SRMX = optimizableVariable('SRMX',[0.1 50],'Type','real');
            
            if Model.Fix_C8
                params = [SD6MM,SD6XX,SD6MX];
            else
                params = [SD6MM,SD6XX,SD6MX,SD8MM,SD8XX,SD8MX];
            end
            
            if ~Model.Fix_Alpha
                params = [params,SAMM,SAXX,SAMX];
            end
            
            params = [params,SRMM,SRXX,SRMX];
        end
        
        if ~Model.Fix_Charge
            params = [params,SQ];
        end
    case {'BH' 'BD' 'BE'}
        if ~isfield(Model,'InnerRange')
            Model.InnerRange = false;
        end
        SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
        
        if Model.SigmaEpsilon
            if Model.InnerRange % For BH model, inner range is when gamma is [0, 6]

                % r_0
                Sr0MM = optimizableVariable('r0_MM',[0 0.15],'Type','real'); % Units: nm
                Sr0XX = optimizableVariable('r0_XX',[0 0.25],'Type','real'); % Units: nm
                Sr0MX = optimizableVariable('r0_MX',[0 0.25],'Type','real'); % Units: nm

                % epsilon
                SepsilonMM = optimizableVariable('epsilon_MM',[50 1e7],'Type','real','Transform','log'); % Units: kJ/mol
                SepsilonXX = optimizableVariable('epsilon_XX',[1e-1 1e6],'Type','real','Transform','log'); % Units: kJ/mol
                SepsilonMX = optimizableVariable('epsilon_MX',[1e-1 1e7],'Type','real','Transform','log'); % Units: kJ/mol

                % gamma
                SgammaMM = optimizableVariable('gamma_MM',[0 3],'Type','real'); % Units: kJ/mol
                SgammaXX = optimizableVariable('gamma_XX',[0 3],'Type','real'); % Units: kJ/mol
                SgammaMX = optimizableVariable('gamma_MX',[0 3],'Type','real'); % Units: kJ/mol
            else
                % r_0
                Sr0MM = optimizableVariable('r0_MM',[0 0.5],'Type','real'); % Units: nm
                Sr0XX = optimizableVariable('r0_XX',[0.15 0.7],'Type','real'); % Units: nm
                Sr0MX = optimizableVariable('r0_MX',[0 0.7],'Type','real'); % Units: nm
                
                % epsilon
                SepsilonMM = optimizableVariable('epsilon_MM',[0 400],'Type','real'); % Units: kJ/mol
                SepsilonXX = optimizableVariable('epsilon_XX',[0 250],'Type','real'); % Units: kJ/mol
                SepsilonMX = optimizableVariable('epsilon_MX',[0 400],'Type','real'); % Units: kJ/mol
                
                % gamma
                SgammaMM = optimizableVariable('gamma_MM',[7 20],'Type','real'); % Units: kJ/mol
                SgammaXX = optimizableVariable('gamma_XX',[7 20],'Type','real'); % Units: kJ/mol
                SgammaMX = optimizableVariable('gamma_MX',[7 20],'Type','real'); % Units: kJ/mol
            end
            if Model.Additivity
                params = [Sr0MM,Sr0XX,SepsilonMM,SepsilonXX,SgammaMX];
            else
                params = [Sr0MM,Sr0XX,Sr0MX,SepsilonMM,SepsilonXX,SepsilonMX,SgammaMM,SgammaXX,SgammaMX];
            end
            
            if ~Model.Fix_Charge
                params = [params,SQ];
            end
        else
            SDMM = optimizableVariable('SDMM',[0 1000],'Type','real');
            SDXX = optimizableVariable('SDXX',[0 40],'Type','real');
            SDMX = optimizableVariable('SDMX',[0 40],'Type','real');

            SRMM = optimizableVariable('SRMM',[0.1 50],'Type','real');
            SRXX = optimizableVariable('SRXX',[0.1 50],'Type','real');
            SRMX = optimizableVariable('SRMX',[0.1 50],'Type','real');

            SDMM2= optimizableVariable('SDMM2',[0 1000],'Type','real');

            % Alpha exponential parameters
            SAMM = optimizableVariable('SAMM',[0.1 10],'Type','real');
            SAXX = optimizableVariable('SAXX',[0.1 10],'Type','real');
            SAMX = optimizableVariable('SAMX',[0.1 10],'Type','real');

            if Model.Additivity
                params = [SDMM,SDXX,SRMM,SRXX];
            else
                params = [SDMM,SDXX,SDMX,SRMM,SRXX,SRMX];
            end

            if ~Model.Fix_Alpha
                params = [params,SAMM,SAXX,SAMX];
            end

            if ~Model.Fix_Charge
                params = [params,SQ];
            end

            if Model.Additivity && Model.Additional_MM_Disp
                params = [params,SDMM2];
            end
        end    
    case 'BF'
        SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
        
        % sigma
        SsigmaMM = optimizableVariable('sigma_MM',[0.05 0.5],'Type','real'); % Units: nm
        SsigmaXX = optimizableVariable('sigma_XX',[0.15 0.7],'Type','real'); % Units: nm
        SsigmaMX = optimizableVariable('sigma_MX',[0.1 0.6],'Type','real'); % Units: nm

        % epsilon
        SepsilonMM = optimizableVariable('epsilon_MM',[0 1000],'Type','real'); % Units: kJ/mol
        SepsilonXX = optimizableVariable('epsilon_XX',[0 220],'Type','real'); % Units: kJ/mol
        SepsilonMX = optimizableVariable('epsilon_MX',[0 400],'Type','real'); % Units: kJ/mol

        % gamma
        SgammaMM = optimizableVariable('gamma_MM',[0 10],'Type','real'); % Units: kJ/mol
        SgammaXX = optimizableVariable('gamma_XX',[0 10],'Type','real'); % Units: kJ/mol
        SgammaMX = optimizableVariable('gamma_MX',[0 10],'Type','real'); % Units: kJ/mol
        
        if Model.Additivity
            params = [SsigmaMM,SsigmaXX,SepsilonMM,SepsilonXX,SgammaMX];
        else
            params = [SsigmaMM,SsigmaXX,SsigmaMX,SepsilonMM,SepsilonXX,SepsilonMX,SgammaMM,SgammaXX,SgammaMX];
        end

        if ~Model.Fix_Charge
            params = [params,SQ];
        end
    case {'JC' 'Mie'}
        if Model.SigmaEpsilon
            SDMM = optimizableVariable('sigma_MM',[0.05 0.46],'Type','real'); % Units: nm
            SDXX = optimizableVariable('sigma_XX',[0.13 0.55],'Type','real'); % Units: nm
            SDMX = optimizableVariable('sigma_MX',[0 0.6],'Type','real'); % Units: nm
            
            SRMM = optimizableVariable('epsilon_MM',[0 300],'Type','real'); % Units: kJ/mol
            SRXX = optimizableVariable('epsilon_XX',[0 50],'Type','real'); % Units: kJ/mol
            SRMX = optimizableVariable('epsilon_MX',[0 100],'Type','real'); % Units: kJ/mol
            
            SDMM2= optimizableVariable('epsilon_MM2',[0 1000],'Type','real'); % Units: kJ/mol
            
        else
            SDMM = optimizableVariable('SDMM',[0 1000],'Type','real');
            SDXX = optimizableVariable('SDXX',[0 40],'Type','real');
            SDMX = optimizableVariable('SDMX',[0 40],'Type','real');
            
            SRMM = optimizableVariable('SRMM',[0.1 50],'Type','real');
            SRXX = optimizableVariable('SRXX',[0.1 50],'Type','real');
            SRMX = optimizableVariable('SRMX',[0.1 50],'Type','real');

            SDMM2= optimizableVariable('SDMM2',[0 1000],'Type','real');
        end
        
        SnMM = optimizableVariable('n_MM',[6 12],'Type','real');
        SnXX = optimizableVariable('n_XX',[6 12],'Type','real');
        SnMX = optimizableVariable('n_MX',[6 12],'Type','real');
        
        if Model.Fix_Charge
            if Model.Additivity
                if Model.Additional_MM_Disp
                    params = [SDMM,SDXX,SRMM,SRXX,SDMM2];
                else
                    params = [SDMM,SDXX,SRMM,SRXX];
                end
            else
                params = [SDMM,SDXX,SDMX,SRMM,SRXX,SRMX];
            end
        else
            SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
            if Model.Additivity
                if Model.Additional_MM_Disp
                    params = [SDMM,SDXX,SRMM,SRXX,SQ,SDMM2];
                else
                    params = [SDMM,SDXX,SRMM,SRXX,SQ];
                end
            else
                params = [SDMM,SDXX,SDMX,SRMM,SRXX,SRMX,SQ];
            end
        end
        if strcmp(Model.Theory,'Mie') && ~Model.Fix_Mie_n
            if Model.Additivity
                params = [params,SnMX];
            else
                params = [params,SnMM,SnXX,SnMX];
            end
        end
    end

    if ~isempty(Model.Additional_GAdjust)
        % GAdjust are N x 3 arrays of gaussian parameters
        % (i , 1) is the Gaussian height of the ith adjustment (may be negative or
        % positive)
        % (i , 2) is the center point of the ith Gaussian (should be positive)
        % (i , 3) is the standard deviation or width (negative and positive values
        % are the same)
        for idx = 1:length(Model.Additional_GAdjust)
            int = [Model.Additional_GAdjust{idx} '_' num2str(idx)];
            GA = optimizableVariable(['GA_' int],Model.Additional_GAdjust_Ranges{idx}(1,:),'Type','real'); % Gaussian depth in kJ/mol (should be negative)
            GB = optimizableVariable(['GB_' int],Model.Additional_GAdjust_Ranges{idx}(2,:),'Type','real'); % Gaussian center in nm (should be positive)
            GC = optimizableVariable(['GC_' int],Model.Additional_GAdjust_Ranges{idx}(3,:),'Type','real'); % Gaussian width in nm(positive)
            params = [params,GA,GB,GC];
        end
    end
    
    if ~isempty(Model.Additional_Function.MM.Range) && Model.Additional_Function.MM.N >= 0
        SNMM = optimizableVariable('SNMM',Model.Additional_Function.MM.Range,'Type','real');
        params = [params,SNMM];
    end
    if ~isempty(Model.Additional_Function.XX.Range) && Model.Additional_Function.XX.N >= 0
        SNXX = optimizableVariable('SNXX',Model.Additional_Function.XX.Range,'Type','real');
        params = [params,SNXX];
    end
    if ~isempty(Model.Additional_Function.MX.Range) && Model.Additional_Function.MX.N >= 0
        SNMX = optimizableVariable('SNMX',Model.Additional_Function.MX.Range,'Type','real');
        params = [params,SNMX];
    end
    
end