% Self contained Theory selection subroutine for CRYSTAL17
function [DFT,disptext,dispswitch,parameters,version] = Select_Theory_for_CRYSTAL17(theory,dispsn,LIMBEK,current_crystal_type,GridSize)

    if strcmp(current_crystal_type,'Atom')
        DFTspin = [newline 'SPIN' newline];
        Spinlock = [newline 'SPINLOCK' newline '1 -7' newline];
        HFspin = 'UHF';
    else
        DFTspin = newline;
        Spinlock = newline;
        HFspin = 'RHF';
    end

    par_input_d3bj = false;
    par_input_d30 = false;
    version = '';
    if strcmp(dispsn,'D3') || strcmp(dispsn,'D3TB')
        dispersion = '-D3';
        disptext = ['_' dispsn];
        dispswitch = true;
    elseif strcmp(dispsn,'D2')
        dispersion = '';
        disptext = '';
        dispswitch = true;
    else
        dispersion = '';
        disptext = '';
        dispswitch = false;
    end

    % Select Theory Options
    if strcmp(theory,'HF')
        DFT = [HFspin Spinlock]; % RHF is the default
                
        S6 = 1.0000;
        S8 = 0.9171;
        A1 = 0.3385;
        A2 = 2.8830;
        par_input_d3bj = true;
        
    elseif strcmp(theory,'B3LYP')
        DFT = ['DFT' DFTspin 'B3LYP' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if strcmp(dispsn,'D3M')
            S6 = 1.0000;
            S8 = 1.4667;
            A1 = 0.2787;
            A2 = 4.6063;
            par_input_d3bj = true;
        end
        
	elseif strcmp(theory,'SVWN') || strcmp(theory,'LSDA')
        DFT = ['DFT' DFTspin 'SVWN' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
	elseif strcmp(theory,'HSE06')
        DFT = ['DFT' DFTspin 'HSE06' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
	elseif strcmp(theory,'HSEsol')
        DFT = ['DFT' DFTspin 'HSEsol' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
	elseif strcmp(theory,'B97')
        DFT = ['DFT' DFTspin 'B97' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
	elseif strcmp(theory,'wB97')
        DFT = ['DFT' DFTspin 'wB97' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
	elseif strcmp(theory,'wB97X')
        DFT = ['DFT' DFTspin 'wB97X' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
	elseif strcmp(theory,'PBE')
        if isempty(dispersion)
            DFT = ['DFT' DFTspin 'PBEXC' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        else
            DFT = ['DFT' DFTspin 'PBE' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        end
        
	elseif strcmp(theory,'LC-wPBE')
        DFT = ['DFT' DFTspin 'LC-wPBE' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
	elseif strcmp(theory,'LC-wPBEsol')
        DFT = ['DFT' DFTspin 'LC-wPBEsol' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
	elseif strcmp(theory,'PBE0')
        DFT = ['DFT' DFTspin 'PBE0' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
    
    elseif strcmp(theory,'PBESOL')
        DFT = ['DFT' DFTspin 'PBESOLXC' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        S6 = 1.000;
        S8 = 2.9491;
        A1 = 0.4466;
        A2 = 6.1742;
        par_input_d3bj = true;
        
    elseif strcmp(theory,'PBESOL0')
        DFT = ['DFT' DFTspin 'PBESOL0' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'BLYP')
        DFT = ['DFT' DFTspin 'BLYP' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
    elseif strcmp(theory,'LC-BLYP')
        DFT = ['DFT' DFTspin 'LC-BLYP' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'LC-wBLYP')
        DFT = ['DFT' DFTspin 'LC-wBLYP' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'CAM-B3LYP')
        DFT = ['DFT' DFTspin 'CAM-B3LYP' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        S6 = 1.000;
        S8 = 2.0674;
        A1 = 0.3708;
        A2 = 5.4743;
        par_input_d3bj = true;
        
    elseif strcmp(theory,'PW1PW')
        DFT = ['DFT' DFTspin 'PW1PW' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
    elseif strcmp(theory,'M06')
        DFT = ['DFT' DFTspin 'M06' dispersion newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        
    elseif strcmp(theory,'RSHXLDA')
        DFT = ['DFT' DFTspin 'RSHXLDA' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'HISS')
        DFT = ['DFT' DFTspin 'HISS' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'SC-BLYP')
        DFT = ['DFT' DFTspin 'SC-BLYP' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'M06L')
        DFT = ['DFT' DFTspin 'M06L' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        S6 = 1.0000;
        S8 = 0.0000;
        rs6 = 1.5810;
        rs8 = 1.0000;
        alpha6 = 14.0000;
        alpha8 = 16.0000;
        par_input_d30 = true;
        
    elseif strcmp(theory,'M05')
        DFT = ['DFT' DFTspin 'M05' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'M052X')
        DFT = ['DFT' DFTspin 'M052X' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end
        
    elseif strcmp(theory,'M062X')
        DFT = ['DFT' DFTspin 'M062X' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        S6 = 1.0000;
        S8 = 0.0000;
        rs6 = 1.6190;
        rs8 = 1.0000;
        alpha6 = 14.0000;
        alpha8 = 16.0000;
        par_input_d30 = true;
        
    elseif strcmp(theory,'M06HF')
        DFT = ['DFT' DFTspin 'M06HF' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'SOGGAXC')
        DFT = ['DFT' DFTspin 'SOGGAXC' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'B3PW')
        DFT = ['DFT' DFTspin 'B3PW' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'B1WC')
        DFT = ['DFT' DFTspin 'B1WC' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'WC1LYP')
        DFT = ['DFT' DFTspin 'WC1LYP' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'B97H')
        DFT = ['DFT' DFTspin 'B97H' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    elseif strcmp(theory,'PBE0-13')
        DFT = ['DFT' DFTspin 'PBE0-13' newline GridSize newline 'LIMBEK' newline LIMBEK newline 'END[DFT]' Spinlock];
        
        if ~isempty(dispersion)
            warning(['Undefined Empirical Dispersion Parameters for theory type: ' ...
                theory '. Excluding dispersion correction from this calculation.'])
            disptext = '';
            dispswitch=false;
        end 
        
    else
        error(['Unknown Theory Type: ' theory])
    end
    
    if par_input_d3bj && dispswitch
        parameters = ['S6' newline num2str(S6) newline 'S8' newline num2str(S8) newline ...
            'A1' newline num2str(A1) newline 'A2' newline num2str(A2) newline 'PRINTC6'];
        version = ['VERSION' newline '4'];
    elseif par_input_d30 && dispswitch
        parameters = ['S6' newline num2str(S6) newline 'S8' newline num2str(S8) newline ...
            'RS6' newline num2str(rs6) newline 'RS8' newline num2str(rs8) newline ...
            'ALPHA6' newline num2str(alpha6) newline 'ALPHA8' newline num2str(alpha8) newline 'PRINTC6'];
        version = ['VERSION' newline '3'];
    else
        parameters = '';
    end
end