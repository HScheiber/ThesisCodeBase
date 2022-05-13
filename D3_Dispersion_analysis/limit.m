function [iat,jat,iadr,jadr] = limit(iat,jat)
      %integer iat,jat,iadr,jadr,i
      iadr=1;
      jadr=1;
      i=100;
      
    while (iat > 100) 
        iat = iat-100;
        iadr = iadr+1;
    end

	while (jat > 100)
        jat=jat-100;
        jadr=jadr+1;
	end

 end