% Copy the C6AB function from memory, with reference C6 data
function [c6ab,maxci] = copyc6
            
    max_elem = 94;
    minc6 = false;
    maxc6 = false;
    maxc = 5;
    minc6list = false;
    maxc6list = false;
    load('pars.mat','nlines','pars')
    maxci = zeros(max_elem,1);
    c6ab = -1.*ones(max_elem,max_elem,maxc,maxc,3);
    
    % read file
    kk=1;
	if (minc6 || maxc6)   %only use values for cn=minimum
        for i=1:94
            if (minc6list(i))
                c6ab(i,:,1,:,2)=10000000.0;
                c6ab(:,i,:,1,3)=10000000.0;
            end
        end


        for nn=1:nlines
            special=false;
            iat=int(pars(kk+1));
            jat=int(pars(kk+2));
            call limit(iat,jat,iadr,jadr)

            if (minc6list(iat))  %only CN=minimum for iat
                special=true;
                maxci(iat)=1;
                maxci(jat)=max(maxci(jat),jadr);

                if(pars(kk+3).le.c6ab(iat,jat,1,jadr,2))

                    c6ab(iat,jat,1,jadr,1)=pars(kk);
                    c6ab(iat,jat,1,jadr,2)=pars(kk+3);
                    c6ab(iat,jat,1,jadr,3)=pars(kk+4);

                    c6ab(jat,iat,jadr,1,1)=pars(kk);
                    c6ab(jat,iat,jadr,1,2)=pars(kk+4);
                    c6ab(jat,iat,jadr,1,3)=pars(kk+3);
                end
            end

            if (minc6list(jat))  %only CN=minimum for jat
                special=true;
                maxci(iat)=max(maxci(iat),iadr);
                maxci(jat)=1;

                if (pars(kk+4).le.c6ab(iat,jat,iadr,1,3))

                    c6ab(iat,jat,iadr,1,1)=pars(kk); 
                    c6ab(iat,jat,iadr,1,2)=pars(kk+3);
                    c6ab(iat,jat,iadr,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,iadr,1)=pars(kk);
                    c6ab(jat,iat,1,iadr,2)=pars(kk+4);
                    c6ab(jat,iat,1,iadr,3)=pars(kk+3);
                end
            end


            if (minc6list(iat).and.minc6list(jat))  %only CN=minimum for iat and jat
                special = true;
                maxci(jat)=1;
                maxci(iat)=1;

                if(pars(kk+4).le.c6ab(iat,jat,1,1,3) && ...
                    pars(kk+3).le.c6ab(iat,jat,1,1,2))

                    c6ab(iat,jat,1,1,1)=pars(kk);  
                    c6ab(iat,jat,1,1,2)=pars(kk+3);
                    c6ab(iat,jat,1,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,1,1)=pars(kk);
                    c6ab(jat,iat,1,1,2)=pars(kk+4);
                    c6ab(jat,iat,1,1,3)=pars(kk+3);
                end
            end



            if (maxc6list(iat)) %only CN=maximum for iat
                special = true;

                maxci(iat)=1;
                maxci(jat)=max(maxci(jat),jadr);

                if(pars(kk+3).ge.c6ab(iat,jat,1,jadr,2))

                    c6ab(iat,jat,1,jadr,1)=pars(kk);  
                    c6ab(iat,jat,1,jadr,2)=pars(kk+3);
                    c6ab(iat,jat,1,jadr,3)=pars(kk+4);

                    c6ab(jat,iat,jadr,1,1)=pars(kk);
                    c6ab(jat,iat,jadr,1,2)=pars(kk+4);
                    c6ab(jat,iat,jadr,1,3)=pars(kk+3);
                end
            end
            if (maxc6list(jat)) %only CN=maximum for jat
                special=true;

                maxci(jat)=1;
                maxci(iat)=max(maxci(iat),iadr);

                if(pars(kk+4).ge.c6ab(iat,jat,iadr,1,3))

                    c6ab(iat,jat,iadr,1,1)=pars(kk);
                    c6ab(iat,jat,iadr,1,2)=pars(kk+3);
                    c6ab(iat,jat,iadr,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,iadr,1)=pars(kk);
                    c6ab(jat,iat,1,iadr,2)=pars(kk+4);
                    c6ab(jat,iat,1,iadr,3)=pars(kk+3);
                end
            end

            if (maxc6list(iat).and.maxc6list(jat))  %only CN=maximum for iat and jat
                special=true;
                maxci(jat)=1;
                maxci(iat)=1;

                if(pars(kk+4).ge.c6ab(iat,jat,1,1,3) && ...
                    pars(kk+3).ge.c6ab(iat,jat,1,1,2))

                    c6ab(iat,jat,1,1,1)=pars(kk);  
                    c6ab(iat,jat,1,1,2)=pars(kk+3);
                    c6ab(iat,jat,1,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,1,1)=pars(kk); 
                    c6ab(jat,iat,1,1,2)=pars(kk+4);
                    c6ab(jat,iat,1,1,3)=pars(kk+3);
                end
            end

            if (minc6list(iat).and.maxc6list(jat))  %only CN=minimum for iat 
                special = true;                            %and CN=maximum jat
                maxci(jat)=1;
                maxci(iat)=1;

                if(pars(kk+4).ge.c6ab(iat,jat,1,1,3) && ...
                    pars(kk+3).le.c6ab(iat,jat,1,1,2))

                    c6ab(iat,jat,1,1,1)=pars(kk);
                    c6ab(iat,jat,1,1,2)=pars(kk+3);
                    c6ab(iat,jat,1,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,1,1)=pars(kk); 
                    c6ab(jat,iat,1,1,2)=pars(kk+4);
                    c6ab(jat,iat,1,1,3)=pars(kk+3);
                end
            end

            if (maxc6list(iat).and.minc6list(jat))  %only CN=maximum for iat
                special = true;                             %  and CN=minimum for jat
                maxci(jat) = 1;
                maxci(iat) = 1;

                if(pars(kk+4).le.c6ab(iat,jat,1,1,3) && ...
                    pars(kk+3).ge.c6ab(iat,jat,1,1,2))

                    c6ab(iat,jat,1,1,1)=pars(kk);
                    c6ab(iat,jat,1,1,2)=pars(kk+3);
                    c6ab(iat,jat,1,1,3)=pars(kk+4);

                    c6ab(jat,iat,1,1,1)=pars(kk);
                    c6ab(jat,iat,1,1,2)=pars(kk+4);
                    c6ab(jat,iat,1,1,3)=pars(kk+3);
                end
            end

            if (~special)

                maxci(iat)=max(maxci(iat),iadr);
                maxci(jat)=max(maxci(jat),jadr);

                c6ab(iat,jat,iadr,jadr,1)=pars(kk);
                c6ab(iat,jat,iadr,jadr,2)=pars(kk+3);
                c6ab(iat,jat,iadr,jadr,3)=pars(kk+4);

                c6ab(jat,iat,jadr,iadr,1)=pars(kk);
                c6ab(jat,iat,jadr,iadr,2)=pars(kk+4);
                c6ab(jat,iat,jadr,iadr,3)=pars(kk+3);
            end
                kk=(nn*5)+1;
        end

    else %no min/max at all 

        for nn = 1:nlines
            iat=pars(kk+1);
            jat=pars(kk+2);
            [iat,jat,iadr,jadr] = limit(iat,jat);
            maxci(iat) = max(maxci(iat),iadr);
            maxci(jat) = max(maxci(jat),jadr);

            c6ab(iat,jat,iadr,jadr,1)=pars(kk);
            c6ab(iat,jat,iadr,jadr,2)=pars(kk+3);
            c6ab(iat,jat,iadr,jadr,3)=pars(kk+4);

            c6ab(jat,iat,jadr,iadr,1)=pars(kk);
            c6ab(jat,iat,jadr,iadr,2)=pars(kk+4);
            c6ab(jat,iat,jadr,iadr,3)=pars(kk+3);
            kk=(nn*5)+1;
        end
	end
end 
