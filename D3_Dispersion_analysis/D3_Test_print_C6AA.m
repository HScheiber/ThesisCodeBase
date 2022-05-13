
echo = true;
n = 4;
iz = [3 3 9 9];

if echo
    ida = zeros(94,1);
    for i=1:n
        ida(iz(i)) = ida(iz(i))+1;
    end
    disp('C6 coefficients used:')
	for i=1:94
    	if(ida(i) > 0)
        	disp([num2str(mxc(i)) ' C6 for element ' num2str(i)])
            for j = 1:maxc
                if (c6ab(i,i,j,j,1) > 0)
                    disp([ 'Z = ' num2str(i) ' CN = ' num2str(c6ab(i,i,j,j,2),'%6.3f') ...
                        ' C6(AA)= ' num2str(c6ab(i,i,j,j,1),'%8.2f') ])
                end
            end
    	end
	end
end