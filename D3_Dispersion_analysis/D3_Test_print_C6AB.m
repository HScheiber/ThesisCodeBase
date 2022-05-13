% Load c6 values
[c6ab,mxc] = copyc6;
echo = true;
dum = []; % empty numerical value
xyz = [0 3.2310 3.47165; 2.83459 1.63655 8.10051; 0 3.27310 0; 2.83459 1.63655 4.62886];
cn = [13.231 13.231 4.126 4.126];
iz = [3 3 9 9];
autoang = 0.52917726;
n=4;
version = 4; % D3(BJ);
fix = [false false false false];
load('r2r4.mat','r2r4')
r0ab = setr0ab;
TZ = false;
func = 'PBE';
fdum = false;
[s6,rs6,s18,rs18,alp] = setfuncpar(func,version,TZ);

if echo
	disp(['#               XYZ [au]  ' repmat(' ', 1, 12) ...
    	'R0(AA) [Ang.]' repmat(' ', 1, 2) 'CN' repmat(' ', 1, 7) ...
    	'C6(AA)     C8(AA)   C10(AA) [au] '])
    x=0;
    btmp='';
    
    for i = 1:n
        z = iz(i);
        c6 = getc6(c6ab,mxc,iz(i),iz(i),cn(i),cn(i));
        for j=1:n
            dum = getc6(c6ab,mxc,iz(i),iz(j),cn(i),cn(j));
            x = x+dum;
        end
        % compute C8/C10 for output
        c8 = r2r4(iz(i))^2*3.0d0*c6;
        c10 = (49.0d0/40.0d0)*c8^2/c6;
        dum = 0.5*autoang*r0ab(z,z);
        if((version == 4) || (version == 6))
            dum = rs6*0.5*autoang*sqrt(c8/c6);
        end
        atmp = ' ';
        if(fix(i)) 
            atmp = 'f';
            btmp = 'f';
        end
        disp([sprintf('%4d',i) sprintf('%10.5f',xyz(i,1:3)) '   ' ...
            sprintf('%2s',elements('atomic_number',z,'Symbol')) ...
            ' ' sprintf('%2s',atmp) sprintf('%7.3f',dum) '  ' sprintf('%7.3f',cn(i)) ...
            sprintf('%12.1f',c6) sprintf('%12.1f',c8) sprintf('%12.1f',c10)]);
    end
    disp(['molecular C6(AA) [au] = ' sprintf('%12.2f',x)])
    if strcmp(btmp,'f')
        disp('  ')
        disp('Caution: Some coordinates fixed in gradient (marked f, see above).')
        disp('  ')
    end
    if fdum
        disp('Caution: Dummy atoms found and ignored.');
    end
end