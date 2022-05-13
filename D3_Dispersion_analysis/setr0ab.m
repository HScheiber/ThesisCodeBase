% Generate bond distances for all 94 x 94 pairs
function r = setr0ab
    load('r0ab.mat','r0ab')
    max_elem = 94;
    autoang = 0.52917726;
    
    r = nan(max_elem,max_elem);
    k=0;
    for i=1:max_elem
        for j=1:i
            k = k+1;
            r(i,j)=r0ab(k)/autoang;
            r(j,i)=r0ab(k)/autoang;
        end
    end
end