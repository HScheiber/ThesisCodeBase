%% clayff_2004_atom.m
% * This function tries to assign all atoms according to the clayff atom types (with modified atom names by MHolmboe), with some modifications for edges...
% * This clayffmod function was modified accoring to Clays and Clay Minerals, Vol. 64, No. 4, 452?471, 2016. STRUCTURE
% * heal_iterations should be a list of numbers 1-7, corresponding to which
% assigment runs you need in order to heal Al, Mg, Si, O, H and so on...
% * The variables distance_factor and rmaxlong are related to the
% neighbour/bond cutoff radius for each atomtype

%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=clayff_2004_atom(atom,Box_dim)
% # atom=clayff_2004_atom(atom,Box_dim,'clayff_2004','spc')
% # atom=clayff_2004_atom(atom,Box_dim,'clayff_2004','spc',heal_iterations)
%
function atom=clayff_2004_atom(atom,Box_dim,varargin)

if nargin >3
    ffname=varargin{1};
    watermodel=varargin{2};
else
    ffname='clayff_2004';
    watermodel='spc/e';
end

if nargin>4
    heal_iterations=varargin{3};
else
    heal_iterations=1; % No healing, will run much faster
end

atom=element_atom(atom);

distance_factor=1.2;
rmaxlong=2.1;

% Initialize some variables
rm_ind=[];
Heal_Al=0;
Heal_Mgo=0;
Heal_Feo=0;
Heal_Mn=0;
Heal_Si=0;
Heal_O=0;
Heal_H=0; % I do not think this is important
Add_H=0; % Add H to get neutral edge
Add_extra_H=0; % Add H to get positive edge As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
% n=7; % with edge structure
for assignment_run=heal_iterations
    
    % Heal broken bonds in the structure in steps, start from the center of the particle
    if assignment_run==1
        Heal_Al=1
    elseif assignment_run==2
        Heal_Mgo=1
    elseif assignment_run==3
        Heal_Si=1
    elseif assignment_run==4
        Heal_O=1
    elseif assignment_run==5
        Heal_H=1
    elseif assignment_run==6
        Add_H=1
    elseif assignment_run==7
        Add_extra_H=0 % As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
    end
    
    All_Neighbours={};
    Bond_index=zeros(4*size(size(atom,2),1),3);
    Angle_index=zeros(4*size(size(atom,2),1),4);
    
    atom=element_atom(atom);
    
    for i=1:size(atom,2)
        if strncmpi(atom(i).element,{'Si'},2);atom(i).type={'st'};
        elseif strncmpi(atom(i).element,{'Alt'},3);atom(i).type={'at'};
        elseif strncmpi(atom(i).element,{'Al'},2);atom(i).type={'ao'};
        elseif strncmpi(atom(i).element,{'Mg'},2);atom(i).type={'mgo'};
        elseif strncmpi(atom(i).element,{'Fe'},2);atom(i).type={'feo'};
        elseif strncmpi(atom(i).element,{'Ow'},2);atom(i).type={'Ow'};
        elseif strncmpi(atom(i).element,{'Hw'},2);atom(i).type={'Hw'};
        elseif strncmpi(atom(i).element,{'O'},1);atom(i).type={'o'};
        elseif strncmpi(atom(i).element,{'H'},1);atom(i).type={'ho'};
        else
            [atom(i).type]=atom(i).element;
        end
    end
    
    [atom.fftype]=atom.type;
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    Radiiproperties=load('Revised_Shannon_radii.mat');
    %     atom = rmfield(atom,'neigh');
    %     atom=bond_valence_atom(atom,Box_dim,1.25,rmax);
    %     clayff_2004_param(sort(unique([atom.type])),'SPC/E');
    
    if size(atom,2)>5000
        dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,rmaxlong);
    else
        dist_matrix = dist_matrix_atom(atom,Box_dim,1.25,2.25);
    end
    
    XYZ_radii=zeros(length(XYZ_labels),1);
    XYZ_formalcharge=zeros(length(XYZ_labels),1);
    Atom_label=sort(unique([atom.element]));
    for i=1:length(Atom_label)
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
        catch
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
        end
        XYZ_radii(ismember([atom.element],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
        XYZ_formalcharge(ismember([atom.element],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
    end
    
    assignin('caller','XYZ_radii',XYZ_radii);
    assignin('caller','XYZ_formalcharge',XYZ_formalcharge);
    
    
    XYZ_radii(XYZ_radii==0)=distance_factor;
    radius_matrix=repmat(XYZ_radii,1,length(XYZ_radii));
    radius_limit=(radius_matrix+radius_matrix')*distance_factor;
    dist_matrix(dist_matrix==0)=1000;
    bond_matrix=dist_matrix-radius_limit;
    dist_matrix(dist_matrix==1000)=0;
    bond_matrix(bond_matrix>0)=0;
    bond_matrix(bond_matrix<0)=1;
    disp('Radii+Radii limits')
    unique(XYZ_labels,'stable')
    unique(XYZ_radii,'stable')
    unique(radius_limit,'stable')
    assignin('caller','radius_limit',radius_limit);
    assignin('caller','bond_matrix',bond_matrix);
    assignin('caller','dist_matrix',dist_matrix);
    
    atom = rmfield(atom,'neigh');
    for i=1:length(XYZ_labels)
        k=1;j=1;
        bond_ind=find(bond_matrix(:,i));
        while j <= numel(bond_ind) && k <= numel(bond_ind) %<= neigh %atom(i).neigh=[];
            if bond_matrix(bond_ind(j),i)==1
                if XYZ_formalcharge(i)*XYZ_formalcharge(bond_ind(j))<0 && [atom(i).molid] == [atom(bond_ind(j)).molid]
                    [atom(i).neigh.dist(k)]=dist_matrix(bond_ind(j),i);
                    [atom(i).neigh.index(k)]=bond_ind(j);
                    [atom(i).neigh.type(k)]=XYZ_labels(bond_ind(j));
                    [atom(i).neigh.coords(k,:)]=[XYZ_data(bond_ind(j),1) XYZ_data(bond_ind(j),2) XYZ_data(bond_ind(j),3)];
                    [atom(i).neigh.r_vec(k,:)]=[X_dist(bond_ind(j),i) Y_dist(bond_ind(j),i) Z_dist(bond_ind(j),i)];
                    k=k+1;
                end
            end
            j=j+1;
        end
    end
    
    assignin('caller','Bond_index',Bond_index);
    
    b=1;a=1;i=1;nH_extra=0;
    while i <= size(atom,2)
        if mod(i,100)==1
            i-1
        end
        
        if strncmpi([atom(i).resname],'SOL',3)==0 && strncmpi([atom(i).resname],'ION',3)==0
            Neigh_ind=zeros(12,1);
            Neigh_vec=zeros(12,3);
            n=1;neigh=1;
            
            if numel([atom(i).neigh])>0
                index=[atom(i).neigh.index];
                if sum(index)>0
                    for k=1:length(index)
                        j=atom(i).neigh.index(k);
                        r=atom(i).neigh.dist(k);
                        
                        if i < size(radius_limit,2)
                            max_distance=radius_limit(j,i);
                        else
                            disp('Manually setting max_distance to...')
                            max_distance=1.8
                        end
                        if r > 0 && r < max_distance/3
                            r
                            disp('Atoms too close')
                            format compact;
                            format short;
                            i
                            j
                            atom(i).type
                            atom(j).type
                            [atom(i).x atom(i).y atom(i).z]
                            [atom(j).x atom(j).y atom(j).z]
                        end
                        
                        if r > max_distance/3 && r < max_distance
                            
                            neigh=neigh+1;
                            
                            if i < j
                                Bond_index(b,1)= i;
                                Bond_index(b,2)= j;
                            else
                                Bond_index(b,1)= j;
                                Bond_index(b,2)= i;
                            end
                            
                            Bond_index(b,3)= r;
                            Neigh_ind(n,1)= j;
                            Neigh_vec(n,1:3)= [atom(i).neigh.r_vec(k,1) atom(i).neigh.r_vec(k,2) atom(i).neigh.r_vec(k,3)];
                            % Neigh_vec(n,1:3)= [atom(i).neigh(k).r_vec(1) atom(i).neigh(k).r_vec(2) atom(i).neigh(k).r_vec(3)];
                            b=b+1;
                            n=n+1;
                        end
                    end
                end
                
                Neigh_cell = sort([atom(i).neigh.type]);
                if length(Neigh_cell) > 0; Neighbours=strcat(Neigh_cell{:});
                else Neighbours={'Nan'}; end
                Neigh_ind(~any(Neigh_ind,2),:) = [];
                Neigh_vec(~any(Neigh_vec,2),:) = [];
                
                size(Neigh_ind,1);
                for v=1:size(Neigh_ind,1)
                    for w=1:size(Neigh_ind,1)
                        angle=rad2deg(atan2(norm(cross(Neigh_vec(v,:),Neigh_vec(w,:))),dot(Neigh_vec(v,:),Neigh_vec(w,:))));
                        if angle > 0 && angle < 160
                            if v < w
                                Angle_index(a,1)= Neigh_ind(v,1);
                                Angle_index(a,2)= i;
                                Angle_index(a,3)= Neigh_ind(w,1);
                                Angle_index(a,4)= angle;
                                a=a+1;
                            else
                                Angle_index(a,1)= Neigh_ind(w,1);
                                Angle_index(a,2)= i;
                                Angle_index(a,3)= Neigh_ind(v,1);
                                Angle_index(a,4)= angle;
                                a=a+1;
                            end
                        end
                    end
                end
                
                Bond_index=sortrows(Bond_index);
                Angle_index=sortrows(Angle_index);
                
                % If the element is Si
                if strncmpi(atom(i).type,{'st'},2) % Si
                    
                    Neigh_cell = sort([atom(i).neigh.type]);
                    Neighbours=strcat(Neigh_cell{:});
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 4 % Si O O O O
                        atom(i).fftype={'st'};
                    elseif length(Neigh_ind) > 4
                        disp('Si atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='st^';
                    elseif length(Neigh_ind) == 3
                        disp('Si atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='st_';
                        if Heal_Si == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'o'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                    end
                    %             else
                    %                 Neighbours
                    %                 atom(i)=[];
                    %                 disp('removed atom...')
                    %                 i
                    %             end
                end
                
                % If the element is Al
                if strncmpi(atom(i).type,{'ao'},2) % Al
                    
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                        atom(i).fftype={'ao'};
                    elseif sum(strncmp({'o'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                        atom(i).fftype={'ae'};
                    elseif sum(strncmp({'o'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                        atom(i).fftype={'at'};
                    elseif length(Neigh_ind) > 6
                        disp('Al atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ao';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('Al atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ae_';
                        if Heal_Al == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'o'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        
                        
                    end
                    %             else
                    %                 Neighbours
                    %                 atom(i)=[];
                    %                 disp('removed atom...')
                    %                 i
                    %             end
                    
                end
                
                % If the element is Mg
                if strncmpi(atom(i).type,{'mg'},2) % Mgo
                    
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 6 % Mgo O O O O O O
                        if sum(strncmpi([atom.type],'ao',2))<sum(strncmpi([atom.type],'mg',2))
                            atom(i).fftype={'mgh'};
                        else
                            atom(i).fftype={'mgo'};
                        end
                    elseif length(Neigh_ind) > 6
                        disp('mgo atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='mgo^';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('mgo atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='mgo_';
                        if Heal_Mgo == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'o',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                    
                end
                
                % If the element is Ca
                if strncmpi(atom(i).type,{'cao'},3) % cao
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 2 % Cao O O
                        atom(i).fftype={'cao'};
                    elseif length(Neigh_ind) > 2
                        disp('cao atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Cao^';
                    elseif length(Neigh_ind) < 2
                        disp('ca atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ca';
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                end
                
                % If the element is H
                if strncmp(atom(i).type,{'ho'},1) % H
                    
                    if size(Neigh_cell,1) == 1
                        atom(i).fftype={'ho'};
                    elseif length(Neigh_ind) > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='h^';
                        atom(i).fftype='h';
                    elseif length(Neigh_ind) < 1
                        disp('H atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='h_';
                        if Heal_H == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'o',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        elseif length(Neigh_ind) < 1
                            Neighbours
                            %                     atom(i)=[];
                            %                     disp('removed atom...')
                            i
                            %pause
                        end
                        
                    end
                    
                end
                
                % If the element is O
                if strncmp(atom(i).type,{'o'},1)
                    
                    All_Neighbours=[All_Neighbours;Neighbours];
                    if strncmp(Neighbours,'aoaoho',5)
                        atom(i).fftype={'oh'};
                    elseif strncmp(Neighbours,'aoaost',5)
                        atom(i).fftype={'op'};
                    elseif strncmp(Neighbours,'hohomgo',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'hohomg',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'aohoho',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'aohost',6)
                        atom(i).fftype={'oahs'}; % Al-OH-Si for acidic edge
                    elseif  strncmp(Neighbours,'aohomg',5) % 'aohmgo'
                        atom(i).fftype={'ohs'};
                    elseif strncmp(Neighbours,'aomgst',5)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoho',4)
                        atom(i).fftype={'oah'}; % Al-O-H or Al-O-Si
                    elseif strncmp(Neighbours,'aomgost',6)
                        atom(i).fftype={'obos'};
                    elseif strncmp(Neighbours,'host',4)
                        atom(i).fftype={'ohst'};
                    elseif strncmp(Neighbours,'aoomg',5)
                        atom(i).fftype={'obss'};
                    elseif strncmp(Neighbours,'aoatho',6)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoaoat',6)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoat',4)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'atst',4)
                        atom(i).fftype={'oat'};
                    elseif strncmp(Neighbours,'aost',4) % Al-O-H or Al-O-Si
                        if sum(ismember(find(strcmp([atom.fftype],'at')),[atom(i).neigh.index])) < 1
                            atom(i).fftype={'oas'};
                        elseif sum(ismember(find(strcmp([atom.fftype],'at')),[atom(i).neigh.index])) > 0
                            atom(i).fftype={'oat'};
                        end
%                     elseif strncmp(Neighbours,'aoao',4)
%                         atom(i).fftype={'o'};
                    elseif strncmp(Neighbours,'mgomgost',7)
                        atom(i).fftype={'odsub'};
                    elseif strncmp(Neighbours,'stst',4)
                        atom(i).fftype={'ob'}; % Old version atom(i).fftype={'O'};
                    elseif strncmp(Neighbours,'hoho',4)
                        atom(i).fftype={'o*'};
                        disp('Water?')
                    elseif length(Neigh_ind) > 2
                        disp('O atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='o^';
                    elseif length(Neigh_ind) == 1 || strncmp(Neighbours,'aoao',4) || strncmpi(Neighbours,'aomgo',5)
                        disp('O atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='o_';
                        if strcmp(Neighbours,'st')
                            atom(i).fftype={'ost'};
                        elseif strcmp(Neighbours,'ao')
                            atom(i).fftype={'oao'};
                        else
                            disp('O atom under coordinated')
                            i
                            Neighbours
                            atom(i).fftype={'o_'};
                        end
                        if Heal_O == 1
                            try
                                atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'ho',1),1));
                            catch
                                atom(size(atom,2)+1)=atom(i);
                            end
                            atom(end).type={'ho'};
                            atom(end).fftype={'ho'};
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    elseif length(Neigh_ind) == 0
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                    if strcmp(atom(i).fftype,'oalh') && Add_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'h',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                    % Special thingy to add Na on the edge
                    %             if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 ;%&& nH_extra > nH_added;
                    %                 disp('Adding Na!!!')
                    %                 atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                    %                 atom(end).index=size(atom,2);
                    %                 atom(end).fftype={'Na'};
                    %                 [atom(end).type]=atom(end).fftype;
                    %                 Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                    %                 NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]+3*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)));
                    %                 [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    %             end
                    %%
                    if strcmp(atom(i).fftype,'oalsi') && Add_extra_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic h to oalsi-->oahs')
                        atom(i).fftype={'oahs'};
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'h',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                end
            end
        end
        i=i+1;
        [atom.type]=atom.fftype;
    end
end

if sum(strncmp([atom.type],'Ow',2))>0
    ind_Ow=find(strncmpi([atom.type],'Ow',2));
    ind_Hw=sort([ind_Ow+1 ind_Ow+2]);
    [atom(ind_Hw).type]=deal({'Hw'});
elseif sum(strncmp([atom.type],'OW',2))>0
    ind_Ow=find(strncmpi([atom.type],'OW',2));
    ind_Hw1=sort(ind_Ow+1);
    ind_Hw2=sort(ind_Ow+2);
    [atom(ind_Hw1).type]=deal({'HW1'});
    [atom(ind_Hw2).type]=deal({'HW2'});
end

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

[Y,I]=sort(Angle_index(:,2));
Angle_index=Angle_index(I,:);
Angle_index = unique(Angle_index,'rows','stable');

Bond_index(~any(Bond_index,2),:) = [];
Angle_index(~any(Angle_index,2),:) = [];

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end
ffname
watermodel

try
    %% If the usual atom types
    atom = charge_atom(atom,Box_dim,ffname,watermodel);
    %% If new atom types
%     atom_alternative = charge_atom(atom,Box_dim,ffname,watermodel,'more');
%     assignin('caller','atom_alternative',atom_alternative);
catch
    disp('Could not set the charge...')
end

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),'SPC/E');
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));
% atom = charge_clayff_atom(atom,Box_dim)

