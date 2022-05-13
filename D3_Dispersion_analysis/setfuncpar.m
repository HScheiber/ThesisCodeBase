function [s6,a1,s8,a2,alp] = setfuncpar(func,version,TZ)
  
%       version: 2 is D2,3 is D3(0),4 is D3(BJ) ,5 is D3M(0),6 is D3M(BJ)
%       func is XC name 
%       TZ is logical for choosing TZ set of parameters with D3(0)

    % Initialize
    func = strrep(lower(func),'-','');
    s6 = nan;
    a1 = nan;
    s8 = nan;
    a2 = nan;
    alp = nan;
    
% double hybrid values revised according to procedure in the GMTKN30 paper
	if version == 6
    	s6  =1.0d0;
        alp =14.0d0;
        % BJ damping with parameters from ...
        switch  func
        	case ('b2plyp')
            	a1 =0.486434;
            	s8 =0.672820;
            	a2=3.656466;
            	s6  =0.640000;
            case ('b3lyp')
            	a1 =0.278672;
            	s8 =1.466677;
            	a2=4.606311;
            case ('b97d')
            	a1 =0.240184;
            	s8 =1.206988;
            	a2=3.864426;
            case ('blyp')
            	a1 =0.448486;
            	s8 =1.875007;
            	a2=3.610679;
            case ('bp')
            	a1 =0.821850;
            	s8 =3.140281;
            	a2=2.728151;
            case ('pbe')
            	a1 =0.012092;
            	s8 =0.358940;
            	a2=5.938951;
            case ('pbe0')
            	a1 =0.007912;
            	s8 =0.528823;
            	a2=6.162326;
            case ('lcwpbe')
                a1 =0.563761;
                s8 =0.906564;
                a2=3.593680;
            otherwise
            	error('functional name unknown')
        end
	elseif version == 5
    	s6  = 1.0d0;
        alp = 14.0d0;
        % zero damping with parameters from ...
        switch func
        	case ('b2plyp')
            	a1 =1.313134;
            	s8 =0.717543;
            	a2=0.016035;
            	s6  =0.640000;
            case ('b3lyp')
            	a1 =1.338153;
            	s8 =1.532981;
            	a2=0.013988;
            case ('b97d')
            	a1 =1.151808;
            	s8 =1.020078;
            	a2=0.035964;
            case ('blyp')
            	a1 =1.279637;
            	s8 =1.841686;
            	a2=0.014370;
            case ('bp')
            	a1 =1.233460;
            	s8 =1.945174;
            	a2=0.000000;
            case ('pbe')
            	a1 =2.340218;
            	s8 =0.000000;
            	a2=0.129434;
            case ('pbe0')
            	a1 =2.077949;
            	s8 =0.000081;
            	a2=0.116755;
            case ('lcwpbe')
            	a1 =1.366361;
            	s8 =1.280619;
            	a2=0.003160;
            otherwise
            	error( 'functional name unknown' )
        end
    % DFTD3 with Becke Johnson finite damping, variant 2 with their radii 
    % SE: Alp is only used in 3body calculations
	elseif version == 4

        s6 = 1.0d0;
        alp = 14.0d0;

        switch func
            case ('bp')
                a1 =0.3946;
                s8 =3.2822;
                a2=4.8516;
            case ('blyp')
                a1 =0.4298;
                s8 =2.6996;
                a2=4.2359;
            case ('revpbe')
                a1 =0.5238;
                s8 =2.3550;
                a2=3.5016;
            case ('rpbe')
                a1 =0.1820;
                s8 =0.8318;
                a2=4.0094;
            case ('b97d')
                a1 =0.5545;
                s8 =2.2609;
                a2=3.2297;
            case ('pbe')
                a1 =0.4289;
                s8 =0.7875;
                a2=4.4407;
            case ('rpw86pbe')
                a1 =0.4613;
                s8 =1.3845;
                a2=4.5062;
            case ('b3lyp')
                a1 =0.3981;
                s8 =1.9889;
                a2=4.4211;
            case ('tpss')
                a1 =0.4535;
                s8 =1.9435;
                a2=4.4752;
            case ('hf')
                a1 =0.3385;
                s8 =0.9171;
                a2=2.8830;
            case ('tpss0')
                a1 =0.3768;
                s8 =1.2576;
                a2=4.5865;
            case ('pbe0')
                a1 =0.4145;
                s8 =1.2177;
                a2=4.8593;
            case ('hse06')
                a1 =0.383;
                s8 =2.310;
                a2=5.685;
            case ('revpbe38')
                a1 =0.4309;
                s8 =1.4760;
                a2=3.9446;
            case ('pw6b95')
                a1 =0.2076;
                s8 =0.7257;
                a2=6.3750;
            case ('b2plyp')
                a1 =0.3065;
                s8 =0.9147;
                a2=5.0570;
                s6=0.64d0;
            case ('dsdblyp')
                a1 =0.0000;
                s8 =0.2130;
                a2=6.0519;
                s6=0.50d0;
            case ('dsdblypfc')
                a1 =0.0009;
                s8 =0.2112;
                a2=5.9807;
                s6=0.50d0;
            case ('bop')
                a1 =0.4870;
                s8 =3.2950;
                a2=3.5043;
            case ('mpwlyp')
                a1 =0.4831;
                s8 =2.0077;
                a2=4.5323;
            case ('olyp')
                a1 =0.5299;
                s8 =2.6205;
                a2=2.8065;
            case ('pbesol')
                a1 =0.4466;
                s8 =2.9491;
                a2=6.1742;
            case ('bpbe')
                a1 =0.4567;
                s8 =4.0728;
                a2=4.3908;
            case ('opbe')
                a1 =0.5512;
                s8 =3.3816;
                a2=2.9444;
            case ('ssb')
                a1 =0.0952;
                s8 =0.1744;
                a2=5.2170;
            case ('revssb')
                a1 =0.4720;
                s8 =0.4389;
                a2=4.0986;
            case ('otpss')
                a1 =0.4634;
                s8 =2.7495;
                a2=4.3153;
            case ('b3pw91')
                a1 =0.4312;
                s8 =2.8524;
                a2=4.4693;
            case ('bhlyp')
                a1 =0.2793;
                s8 =1.0354;
                a2=4.9615;
            case ('revpbe0')
                a1 =0.4679;
                s8 =1.7588;
                a2=3.7619;
            case ('tpssh')
                a1 =0.4529;
                s8 =2.2382;
                a2=4.6550;
            case ('mpw1b95')
                a1 =0.1955;
                s8 =1.0508;
                a2=6.4177;
            case ('pwb6k')
                a1 =0.1805;
                s8 =0.9383;
                a2=7.7627;
            case ('b1b95')
                a1 =0.2092;
                s8 =1.4507;
                a2=5.5545;
            case ('bmk')
                a1 =0.1940;
                s8 =2.0860;
                a2=5.9197;
            case ('camb3lyp')
                a1 =0.3708;
                s8 =2.0674;
                a2=5.4743;
            case ('lcwpbe')
                a1 =0.3919;
                s8 =1.8541;
                a2=5.0897;
            case ('b2gpplyp')
                a1 =0.0000;
                s8 =0.2597;
                a2=6.3332;
                s6=0.560;
            case ('ptpss')
                a1 =0.0000;
                s8 =0.2804;
                a2=6.5745;
                s6=0.750;
            case ('pwpb95')
                a1 =0.0000;
                s8 =0.2904;
                a2=7.3141;
                s6=0.820;
            % special HF/DFT with eBSSE correction
            case ('hf/mixed')
                a1 =0.5607;  
                s8 =3.9027;  
                a2=4.5622;  
            case ('hf/sv')
                a1 =0.4249;
                s8 =2.1849;
                a2=4.2783;
            case ('hf/minis')
                a1 =0.1702;
                s8 =0.9841;
                a2=3.8506;
            case ('b3lyp/631gd')
                a1 =0.5014;
                s8 =4.0672;
                a2=4.8409;
            case ('hcth120')
                a1=0.3563;
                s8=1.0821;
                a2=4.3359;
            % DFTB3 old, deprecated parameters:
            %         case ('dftb3')
            %              rs6=0.7461
            %              s18=3.209 
            %              rs18=4.1906
            % special SCCDFTB parametrization
            % full third order DFTB, self consistent charges, hydrogen pair damping with 
            %         exponent 4.2
            case('dftb3')
                a1=0.5719d0;
                s8=0.5883d0;
                a2=3.6017d0;
            case ('pw1pw')
                a1 =0.3807d0;
                s8 =2.3363d0;
                a2=5.8844d0;
            case ('pwgga')
                a1 =0.2211d0;
                s8 =2.6910d0;
                a2=6.7278d0;
            case ('hsesol')
                a1 =0.4650d0;
                s8 =2.9215d0;
                a2=6.2003d0;
            % special HFD3gCPSRB/MINIX parametrization
            case ('hf3c')
                a1=0.4171d0;
                s8=0.8777d0;
                a2=2.9149d0;
            % special HFD3gCPSRB2/ECP2G parametrization
            case ('hf3cv')
                a1=0.3063d0;
                s8=0.5022d0;
                a2=3.9856d0;
            % special PBEhD3gCP/def2mSVP parametrization
            case ('pbeh3c')
                a1=0.4860d0;
                s8=0.0000d0;
                a2=4.5000d0;
            otherwise 
            error('functional name unknown')
        end
    elseif (version == 3)
        % DFTD3
        s6  =1.0d0;
        alp =14.0d0;
        a2=1.0d0;
        % default def2QZVP (almost basis set limit)
        if (~TZ)
            switch func
            	case ('slaterdiracexchange')
            	     a1 =0.999;
            	     s8 =1.957;
            	     a2=0.697;
            	case ('blyp')
            	     a1=1.094;
            	     s8=1.682;
            	case ('bp')
            	     a1=1.139;
            	     s8=1.683;
            	case ('b97d')
            	     a1=0.892;
            	     s8=0.909;
            	case ('revpbe')
            	     a1=0.923;
            	     s8=1.010;
            	case ('pbe')
            	     a1=1.217;
            	     s8=0.722;
            	case ('pbesol')
            	     a1=1.345;
            	     s8=0.612;
            	case ('rpw86pbe')
            	     a1=1.224;
            	     s8=0.901;
            	case ('rpbe')
            	     a1=0.872;
            	     s8=0.514;
            	case ('tpss')
            	     a1=1.166;
            	     s8=1.105;
            	case ('b3lyp')
            	     a1=1.261;
            	     s8=1.703;
            	case ('pbe0')
            	     a1=1.287;
            	     s8=0.928;
            	case ('hse06')
            	     a1=1.129;
            	     s8=0.109;
            	case ('revpbe38')
            	     a1=1.021;
            	     s8=0.862;
            	case ('pw6b95')
            	     a1=1.532;
            	     s8=0.862;
            	case ('tpss0')
            	     a1=1.252;
            	     s8=1.242;
            	case ('b2plyp')
            	     a1=1.427;
            	     s8=1.022;
            	     s6=0.64;
            	case ('pwpb95')
            	     a1=1.557;
            	     s8=0.705;
            	     s6=0.82;
            	case ('b2gpplyp')
            	     a1=1.586;
            	     s8=0.760;
            	     s6=0.56;
            	case ('ptpss')
            	     a1=1.541;
            	     s8=0.879;
            	     s6=0.75;
            	case ('hf')
            	     a1=1.158;
            	     s8=1.746;
            	case ('mpwlyp')
            	     a1=1.239;
            	     s8=1.098;
            	case ('bpbe')
            	     a1=1.087;
            	     s8=2.033;
            	case ('bhlyp')
            	     a1=1.370;
            	     s8=1.442;
            	case ('tpssh')
            	     a1=1.223;
            	     s8=1.219;
            	case ('pwb6k')
            	     a1=1.660;
            	     s8=0.550;
            	case ('b1b95')
            	     a1=1.613;
            	     s8=1.868;
            	case ('bop')
            	     a1=0.929;
            	     s8=1.975;
            	case ('olyp')
            	     a1=0.806;
            	     s8=1.764;
            	case ('opbe')
            	     a1=0.837;
            	     s8=2.055;
            	case ('ssb')
            	     a1=1.215;
            	     s8=0.663;
            	case ('revssb')
            	     a1=1.221;
            	     s8=0.560;
            	case ('otpss')
            	     a1=1.128;
            	     s8=1.494;
            	case ('b3pw91')
            	     a1=1.176;
            	     s8=1.775;
            	case ('revpbe0')
            	     a1=0.949;
            	     s8=0.792;
            	case ('pbe38')
            	     a1=1.333;
            	     s8=0.998;
            	case ('mpw1b95')
            	     a1=1.605;
            	     s8=1.118;
            	case ('mpwb1k')
            	     a1=1.671;
            	     s8=1.061;
            	case ('bmk')
            	     a1=1.931;
            	     s8=2.168;
            	case ('camb3lyp')
            	     a1=1.378;
            	     s8=1.217;
            	case ('lcwpbe')
            	     a1=1.355;
            	     s8=1.279;
            	case ('m05')
            	     a1=1.373;
            	     s8=0.595;
            	case ('m052x')
            	     a1=1.417;
            	     s8=0.000;
            	case ('m06l')
            	     a1=1.581;
            	     s8=0.000;
            	case ('m06')
            	     a1=1.325;
            	     s8=0.000;
            	case ('m062x')
            	     a1=1.619;
            	     s8=0.000;
            	case ('m06hf')
            	     a1=1.446;
            	     s8=0.000;
        % DFTB3 (zeta=4.0), old deprecated parameters
        %         case ('dftb3')
        %              rs6=1.235
        %              s18=0.673
            	case ('hcth120')
            	     a1=1.221;
            	     s8=1.206;
                otherwise
            	     error( 'functional name unknown' )
            end 
        else
            % special TZVPP parameter
            switch func
            	case ('blyp')
            	     a1=1.243;
            	     s8=2.022;
            	case ('bp')
            	     a1=1.221;
            	     s8=1.838;
            	case ('b97d')
            	     a1=0.921;
            	     s8=0.894;
            	case ('revpbe')
            	     a1=0.953;
            	     s8=0.989;
            	case ('pbe')
            	     a1=1.277;
            	     s8=0.777;
            	case ('tpss')
            	     a1=1.213;
            	     s8=1.176;
            	case ('b3lyp')
            	     a1=1.314;
            	     s8=1.706;
            	case ('pbe0')
            	     a1=1.328;
            	     s8=0.926;
            	case ('pw6b95')
            	     a1=1.562;
            	     s8=0.821;
            	case ('tpss0')
            	     a1=1.282;
            	     s8=1.250;
            	case ('b2plyp')
            	     a1=1.551;
            	     s8=1.109;
            	     s6=0.5;
            	otherwise
            	     error( 'functional name unknown (TZ case)' )
            end 
        end
    elseif (version == 2)
        % DFTD2
        a1=1.1d0;
        s8=0.0d0;
        alp=20.0d0;
        switch (func)
            case ('blyp')
                s6=1.2; 
            case ('bp')
            	s6=1.05; 
            case ('b97d')
            	s6=1.25;
            case ('revpbe')
            	s6=1.25;
            case ('pbe')
            	s6=0.75;
            case ('tpss')
            	s6=1.0;
            case ('b3lyp')
            	s6=1.05; 
            case ('pbe0')
            	s6=0.6;  
            case ('pw6b95')
            	s6=0.5;  
            case ('tpss0')
            	s6=0.85;
            case ('b2plyp')
            	s6=0.55;
            case ('b2gpplyp')
            	s6=0.4;
            case ('dsdblyp')
            	s6=0.41;
            	alp=60.0d0;
            otherwise
            	error( 'functional name unknown' )
        end
	end
end