function Text_out = AddCartesianCoord(Template_text,Crystal,idx,FindParam,CoordType)
    if FindParam
        XYZ_Met_Coords = Crystal.FC_Metal{idx};
        XYZ_Hal_Coords = Crystal.FC_Halide{idx};
    else
        XYZ_Met_Coords = Crystal.FC_Metal;
        XYZ_Hal_Coords = Crystal.FC_Halide;
    end

    % Convert lattice param into nm
    a = Crystal.a(idx)/10;
    b = Crystal.b(idx)/10;
    c = Crystal.c(idx)/10;

    % Lattice Parameter matrix
    Lattice_Matrix = [a 0 0; 0 b 0; 0 0 c];

    % Transform LP matrix
    Lattice_Mat_Cartesian = Lattice_Matrix*Crystal.Transform;

    % For gro files
    if strcmp(CoordType,'gro')
        for jj=1:size(XYZ_Met_Coords,1)
            XYZ_Met_Coords(jj,:) = round( ( (XYZ_Met_Coords(jj,:).*[a b c])*Crystal.Transform )*10000)/10000;
            XYZ_Hal_Coords(jj,:) = round( ( (XYZ_Hal_Coords(jj,:).*[a b c])*Crystal.Transform )*10000)/10000;

            Template_text = strrep(Template_text,['##MX' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,1),'%5.3f'));
            Template_text = strrep(Template_text,['##MY' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,2),'%5.3f'));
            Template_text = strrep(Template_text,['##MZ' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,3),'%5.3f'));

            Template_text = strrep(Template_text,['##HX' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,1),'%5.3f'));
            Template_text = strrep(Template_text,['##HY' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,2),'%5.3f'));
            Template_text = strrep(Template_text,['##HZ' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,3),'%5.3f'));
        end

        % Add in lattice Parameters
        Lattice_text = cell(3,3);
        for ii=1:3
            for jj=1:3
                if Lattice_Mat_Cartesian(ii,jj) < 0
                    Lattice_text{ii,jj} = num2str(Lattice_Mat_Cartesian(ii,jj),'%.5f');
                else
                    Lattice_text{ii,jj} = [' ' num2str(Lattice_Mat_Cartesian(ii,jj),'%.5f')];
                end

            end
        end
        
    % For g96 files (extra precision)
    elseif strcmp(CoordType,'g96')
        
        for jj=1:size(XYZ_Met_Coords,1)
            XYZ_Met_Coords(jj,:) = round( ( (XYZ_Met_Coords(jj,:).*[a b c])*Crystal.Transform )*1e9)/1e9;
            XYZ_Hal_Coords(jj,:) = round( ( (XYZ_Hal_Coords(jj,:).*[a b c])*Crystal.Transform )*1e9)/1e9;

            Template_text = strrep(Template_text,['##MX' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,1),'%15.9f'));
            Template_text = strrep(Template_text,['##MY' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,2),'%15.9f'));
            Template_text = strrep(Template_text,['##MZ' num2str(jj) '##'],num2str(XYZ_Met_Coords(jj,3),'%15.9f'));

            Template_text = strrep(Template_text,['##HX' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,1),'%15.9f'));
            Template_text = strrep(Template_text,['##HY' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,2),'%15.9f'));
            Template_text = strrep(Template_text,['##HZ' num2str(jj) '##'],num2str(XYZ_Hal_Coords(jj,3),'%15.9f'));
        end
        
        % Add in lattice Parameters
        Lattice_text = cell(3,3);
        for ii=1:3
            for jj=1:3
                if Lattice_Mat_Cartesian(ii,jj) < 0
                    Lattice_text{ii,jj} = num2str(Lattice_Mat_Cartesian(ii,jj),'%15.9f');
                else
                    Lattice_text{ii,jj} = [' ' num2str(Lattice_Mat_Cartesian(ii,jj),'%15.9f')];
                end

            end
        end  
    end

%            3.000000000000000                   0                   0
%           -1.500000000000000   2.598076211353316                   0
%                    0                           0   4.898979485566356

%    0.30000   0.25981   0.48990   0.00000   0.00000  -0.15000   0.00000  -0.00000  -0.00000
%    v1(x)     v2(y)     v3(z)     v1(y)     v1(z)    v2(x)      v2(z)    v3(x)     v3(y) 

    for i = 1:3
        id = num2str(i);
        Template_text = strrep(Template_text,['##V' id 'X##'],Lattice_text{i,1});
        Template_text = strrep(Template_text,['##V' id 'Y##'],Lattice_text{i,2});
        Template_text = strrep(Template_text,['##V' id 'Z##'],Lattice_text{i,3});
    end
    
    % Output
    Text_out = Template_text;
end