% Function to interpolate C6
function c6 = getc6(c6ab_ref,mxc,ZA,ZB,nci,ncj)

% mxc is an array ordered from Z=1 to Z=94 with how many different reference C6 for each element
%% c6ab_ref is the reference c6 function, a 5D (94 x 94 x 5 x 5 x 3) array
% Data is referenced as follows:
% c6ab_ref(Z_A,Z_B,iA,jB,1) = Stored C6 value between A and B, using the iA
% reference state and the jB reference state.
%% c6ab_ref(Z_A,Z_B,iA,jB,2) = Stored reference Fractional Coordination number of A
%% c6ab_ref(Z_A,Z_B,iA,jB,3) = Stored reference Fractional Coordination number of B

% Ad hock params
    k1 = 16;
    k2 = 4/3;
    k3 = -4;

    % Initialize stuff
    c6mem = -1.d+99;
    rsum = 0.0d0;
    csum = 0.0d0;
    c6 = 0.0d0;
    r_save = 1.0d99;
    
	for i=1:mxc(ZA)
        for j=1:mxc(ZB)
            c6 = c6ab_ref(ZA,ZB,i,j,1);
            
            if(c6 > 0)
            	% c6mem=c6
                cn1 = c6ab_ref(ZA,ZB,i,j,2);
                cn2 = c6ab_ref(ZA,ZB,i,j,3);
                % distance
                r = (cn1-nci)^2 + (cn2-ncj)^2;
                if (r < r_save)
                    r_save = r;
                    c6mem = c6;
                end
                tmp1=exp(k3*r);
                rsum = rsum+tmp1;
                csum = csum+tmp1*c6;
            end
        end
	end

    if(rsum > 1.0d-99)
        c6 = csum/rsum;
    else
        c6 = c6mem;
    end
end 