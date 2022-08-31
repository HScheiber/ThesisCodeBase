function params = bayesopt_params(Model)

    if strcmp(Model.Theory,'TF')
        SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
        
        if Model.SigmaEpsilon
            % r_0
            Sr0MM = optimizableVariable('r0_MM',Model.Sr0MM_Range,'Type','real'); % Units: nm
            Sr0XX = optimizableVariable('r0_XX',Model.Sr0XX_Range,'Type','real'); % Units: nm
            Sr0MX = optimizableVariable('r0_MX',Model.Sr0MX_Range,'Type','real'); % Units: nm
            
            % epsilon
            SepsilonMM = optimizableVariable('epsilon_MM',Model.SepsMM_Range,'Type','real'); % Units: kJ/mol
            SepsilonXX = optimizableVariable('epsilon_XX',Model.SepsXX_Range,'Type','real'); % Units: kJ/mol
            SepsilonMX = optimizableVariable('epsilon_MX',Model.SepsMX_Range,'Type','real'); % Units: kJ/mol
            
            % gamma
            SgammaMM = optimizableVariable('gamma_MM',Model.SgamMM_Range,'Type','real'); % Units: kJ/mol
            SgammaXX = optimizableVariable('gamma_XX',Model.SgamXX_Range,'Type','real'); % Units: kJ/mol
            SgammaMX = optimizableVariable('gamma_MX',Model.SgamMX_Range,'Type','real'); % Units: kJ/mol
            
            if Model.Additivity
                params = [Sr0MM,Sr0XX,SepsilonMM,SepsilonXX,SgammaMX];
            else
                params = [Sr0MM,Sr0XX,Sr0MX,SepsilonMM,SepsilonXX,SepsilonMX,SgammaMM,SgammaXX,SgammaMX];
            end
        else
            % C6 dispersion parameters
            SD6MM = optimizableVariable('SD6MM',Model.SD6MM_Range,'Type','real');
            SD6XX = optimizableVariable('SD6XX',Model.SD6XX_Range,'Type','real');
            SD6MX = optimizableVariable('SD6MX',Model.SD6MX_Range,'Type','real');
            
            % C8 dispersion parameters
            SD8MM = optimizableVariable('SD8MM',Model.SD8MM_Range,'Type','real');
            SD8XX = optimizableVariable('SD8XX',Model.SD8XX_Range,'Type','real');
            SD8MX = optimizableVariable('SD8MX',Model.SD8MX_Range,'Type','real');
            
            % Alpha exponential parameters
            SAMM = optimizableVariable('SAMM',Model.SAMM_Range,'Type','real');
            SAXX = optimizableVariable('SAXX',Model.SAXX_Range,'Type','real');
            SAMX = optimizableVariable('SAMX',Model.SAMX_Range,'Type','real');
            
            % Repulsive wall prefactors
            SRMM = optimizableVariable('SRMM',Model.SRTFMM_Range,'Type','real');
            SRXX = optimizableVariable('SRXX',Model.SRTFXX_Range,'Type','real');
            SRMX = optimizableVariable('SRMX',Model.SRTFMX_Range,'Type','real');
            
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
        
    elseif strcmp(Model.Theory,'BH')
        SQ = optimizableVariable('SQ',Model.Q_Range,'Type','real');
        
        if Model.SigmaEpsilon
            % r_0
            Sr0MM = optimizableVariable('r0_MM',Model.Sr0MM_Range,'Type','real'); % Units: nm
            Sr0XX = optimizableVariable('r0_XX',Model.Sr0XX_Range,'Type','real'); % Units: nm
            Sr0MX = optimizableVariable('r0_MX',Model.Sr0MX_Range,'Type','real'); % Units: nm
            
            % epsilon
            SepsilonMM = optimizableVariable('epsilon_MM',Model.SepsMM_Range,'Type','real'); % Units: kJ/mol
            SepsilonXX = optimizableVariable('epsilon_XX',Model.SepsXX_Range,'Type','real'); % Units: kJ/mol
            SepsilonMX = optimizableVariable('epsilon_MX',Model.SepsMX_Range,'Type','real'); % Units: kJ/mol
            
            % gamma
            SgammaMM = optimizableVariable('gamma_MM',Model.SgamMM_RangeBH,'Type','real'); % Units: kJ/mol
            SgammaXX = optimizableVariable('gamma_XX',Model.SgamXX_RangeBH,'Type','real'); % Units: kJ/mol
            SgammaMX = optimizableVariable('gamma_MX',Model.SgamMX_RangeBH,'Type','real'); % Units: kJ/mol
            
            if Model.Additivity
                params = [Sr0MM,Sr0XX,SepsilonMM,SepsilonXX,SgammaMX];
            else
                params = [Sr0MM,Sr0XX,Sr0MX,SepsilonMM,SepsilonXX,SepsilonMX,SgammaMM,SgammaXX,SgammaMX];
            end
            
            if ~Model.Fix_Charge
                params = [params,SQ];
            end
        else
            SDMM = optimizableVariable('SDMM',Model.SDMM_Range,'Type','real');
            SDXX = optimizableVariable('SDXX',Model.SDXX_Range,'Type','real');
            SDMX = optimizableVariable('SDMX',Model.SDMX_Range,'Type','real');

            SRMM = optimizableVariable('SRMM',Model.SRMM_Range,'Type','real');
            SRXX = optimizableVariable('SRXX',Model.SRXX_Range,'Type','real');
            SRMX = optimizableVariable('SRMX',Model.SRMX_Range,'Type','real');

            SDMM2= optimizableVariable('SDMM2',Model.SDMM2_Range,'Type','real');

            % Alpha exponential parameters
            SAMM = optimizableVariable('SAMM',Model.SAMM_Range,'Type','real');
            SAXX = optimizableVariable('SAXX',Model.SAXX_Range,'Type','real');
            SAMX = optimizableVariable('SAMX',Model.SAMX_Range,'Type','real');

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
        
    else % JC models
        if Model.SigmaEpsilon
            SDMM = optimizableVariable('Sigma_MM',[0.1,0.4],'Type','real'); % Units: nm
            SDXX = optimizableVariable('Sigma_XX',[0.1,0.8],'Type','real'); % Units: nm
            SDMX = optimizableVariable('Sigma_MX',[0.1,0.8],'Type','real'); % Units: nm
            
            if contained_in_cell('MM',Model.Additional_GAdjust) || (Model.Additivity && Model.Additional_MM_Disp)
                SRMM = optimizableVariable('Epsilon_MM',[0,10],'Type','real'); % Units: kJ/mol
            else
                SRMM = optimizableVariable('Epsilon_MM',[0,1000],'Type','real'); % Units: kJ/mol
            end
            SRXX = optimizableVariable('Epsilon_XX',[0,1000],'Type','real'); % Units: kJ/mol
            SRMX = optimizableVariable('Epsilon_MX',[0,1000],'Type','real'); % Units: kJ/mol
            
            SDMM2= optimizableVariable('Epsilon_MM2',[0,1e6],'Type','real'); % Units: kJ/mol
            
        else
            SDMM = optimizableVariable('SDMM',Model.SDMM_Range,'Type','real');
            SDXX = optimizableVariable('SDXX',Model.SDXX_Range,'Type','real');
            SDMX = optimizableVariable('SDMX',Model.SDMX_Range,'Type','real');
            
            SRMM = optimizableVariable('SRMM',Model.SRMM_Range,'Type','real');
            SRXX = optimizableVariable('SRXX',Model.SRXX_Range,'Type','real');
            SRMX = optimizableVariable('SRMX',Model.SRMX_Range,'Type','real');

            SDMM2= optimizableVariable('SDMM2',Model.SDMM2_Range,'Type','real');
        end
        
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