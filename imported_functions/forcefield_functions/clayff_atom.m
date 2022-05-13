%% clayff_atom.m
% * This function tries to assign all atoms according to the clayff atom
% types (with modified atom names by MHolmboe), with some modifications for
% edges...
% * heal_iterations should be a list of numbers 1-7, corresponding to which
% assigment runs you need in order to heal Al, Mg, Si, O, H and so on...
% * The variables distance_factor and rmaxlong are related to the
% neighbour/bond cutoff radius for each atomtype
%
% * This clayff_atom function was modified accoring to Clays and Clay
% Minerals, Vol. 64, No. 4, 452?471, 2016.
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=clayff_atom(atom,Box_dim)
% # atom=clayff_atom(atom,Box_dim,'clayff','spc')
% # atom=clayff_atom(atom,Box_dim,'clayff','spc',heal_iterations)

function atom=clayff_atom(atom,Box_dim,varargin)
%%
format compact;

if nargin >3
    ffname=varargin{1};
    watermodel=varargin{2};
else
    ffname='clayff';
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
        Heal_Al=0
    elseif assignment_run==2
        Heal_Mgo=0
    elseif assignment_run==3
        Heal_Feo=0
    elseif assignment_run==4
        Heal_Mn=0
    elseif assignment_run==5
        Heal_Si=1
    elseif assignment_run==6
        Heal_O=1
    elseif assignment_run==7
        Heal_H=1
    elseif assignment_run==8
        Add_H=1
    elseif assignment_run==9
        Add_extra_H=0 % As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
    end
    
    All_Neighbours={};
    Bond_index=zeros(4*size(size(atom,2),1),3);
    Angle_index=zeros(4*size(size(atom,2),1),4);
    
    [atom.element]=atom.type;
    
    for i=1:size(atom,2)
        if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
        elseif strncmp(atom(i).type,{'SC'},2);atom(i).element={'Si'};
        elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
        elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
        elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
            %         elseif strncmpi(atom(i).type,{'Mn'},2);atom(i).element={'Mn'};
        elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={'Ow'};
        elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={'Hw'};
        elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};
        elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};
        else
            [atom(i).element]=atom(i).type;
        end
    end
    
    [atom.type]=atom.element;
    [atom.fftype]=atom.element;
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    Radiiproperties=load('Revised_Shannon_radii.mat');
%     atom = rmfield(atom,'neigh');
%     atom=bond_valence_atom(atom,Box_dim,1.25,rmax);
%     clayff_param(sort(unique([atom.type])),'SPC/E');
    
    if size(atom,2)>5000
        dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,rmaxlong);
    else
        dist_matrix = dist_matrix_atom(atom,Box_dim,1.25,rmaxlong);
    end
    
    XYZ_radii=zeros(length(XYZ_labels),1);
    XYZ_formalcharge=zeros(length(XYZ_labels),1);
    Atom_label=sort(unique([atom.type]));
    for i=1:length(Atom_label)
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
        catch
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
        end
        XYZ_radii(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
        XYZ_formalcharge(ismember([atom.type],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
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
    %     assignin('caller','newatom',atom);
    
    b=1;a=1;i=1;nH_extra=0;
    while i <= size(atom,2)
        if mod(i,100)==1
            i-1
        end
        
        %         if strncmpi([atom(i).resname],'SOL',3)==0 && strncmpi([atom(i).resname],'ION',3)==0
        
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
                        
                        Neigh_ind(n,1)= j;
                        Neigh_vec(n,1:3)= [atom(i).neigh.r_vec(k,1) atom(i).neigh.r_vec(k,2) atom(i).neigh.r_vec(k,3)];
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
            
            % If the element is Si
            if strncmpi(atom(i).type,{'Si'},2) % Si
                
                Neigh_cell = sort([atom(i).neigh.type]);
                Neighbours=strcat(Neigh_cell{:});
                if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Si O O O O
                    atom(i).fftype={'Si'};
                elseif length(Neigh_ind) > 4
                    disp('Si atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Si^'};
                    atom(i).neigh.dist
                    atom(i).neigh.index
                elseif length(Neigh_ind) < 4 % == 3
                    disp('Si atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Si_'};
                    if Heal_Si == 1
                        atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'O'),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.85*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                else
                    i
                    atom(i).type
                    Neighbours
%                     atom(i)=[];
%                     disp('removed atom...')
%                     i
                end
                %             else
                %                 Neighbours
                %                 atom(i)=[];
                %                 disp('removed atom...')
                %                 i
                %             end
            end
            
            % If the element is Al
            if strncmpi(atom(i).type,{'Al'},2) % Al
                
                if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                    atom(i).fftype={'Al'};
                elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                    atom(i).fftype={'Ale'};
                    %                 elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                    %                     atom(i).fftype={'Alt'};
                elseif length(Neigh_ind) > 6
                    disp('Al atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Al'};
                elseif length(Neigh_ind) == 4
%                     disp('Al atom is tetrahedrally coordinated')
%                     i;
%                     Neighbours;
                    atom(i).fftype={'Alt'};
                elseif length(Neigh_ind) >= 4 && length(Neigh_ind) < 6
                    disp('Al atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Al_'};
                    if Heal_Al == 1
                        disp('Healing Al with O!!!')
                        atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                        atom(end).index=size(atom,2);
                        atom(end).type={'O'};
                        atom(end).fftype={'O'};
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.85*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
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
            if strncmpi(atom(i).type,{'Mg'},2) % Mgo
                
                if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Mgo O O O O O O
                    if sum(strncmpi([atom.type],'Al',2))<sum(strncmpi([atom.type],'Mg',2))
                        atom(i).fftype={'Mgh'};
                    else
                        atom(i).fftype={'Mgo'};
                    end
                elseif length(Neigh_ind) > 6
                    disp('Mgo atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Mgo^'};
                elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                    disp('Mgo atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Mgo_'};
                    if Heal_Mgo == 1
                        atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                        atom(end).index=size(atom,2);
                        atom(end).type={'O'};
                        atom(end).fftype={'O'};
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.7*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
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
            
            % If the element is Mn
            if strncmpi(atom(i).type,{'Mn'},2) % Fe
                
                if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Mn O O O O O O
                    atom(i).fftype={'Mn4'};
                elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Mn O O O O
                    atom(i).fftype={'Mne'};
                elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Mn O O O O
                    atom(i).fftype={'Mnt'};
                elseif length(Neigh_ind) > 6
                    disp('Mn atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Mn^'};
                elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                    disp('Mn atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Mn_'};
                    if Heal_Fe == 1
                        atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                        atom(end).index=size(atom,2);
                        atom(end).type={'O'};
                        atom(end).fftype={'O'};
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.9*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                else
                    Neighbours
                    atom(i).fftype
                    atom(i).type
                    i
                end
                %             else
                %                 Neighbours
                %                 atom(i)=[];
                %                 disp('removed atom...')
                %                 i
                %             end
            end
            
            % If the element is Fe
            if strncmpi(atom(i).type,{'Fe'},2) % Fe
                
                if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Fe O O O O O O
                    atom(i).fftype={'Feo'};
                elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Fe O O O O
                    atom(i).fftype={'Fee'};
                elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Fe O O O O
                    atom(i).fftype={'Fet'};
                elseif length(Neigh_ind) > 6
                    disp('Fe atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Fe^'};
                elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                    disp('Fe atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'Fe_'};
                    if Heal_Fe == 1
                        atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                        atom(end).index=size(atom,2);
                        atom(end).type={'O'};
                        atom(end).fftype={'O'};
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.9*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
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
            
            % If the element is H
            if strncmp(atom(i).type,{'H'},1) % H
                
                if size(Neigh_cell,1) == 1
                        atom(i).fftype={'H'};
                elseif length(Neigh_ind) > 1
                    disp('H atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'H^'};
                    atom(i).fftype={'H'};
                elseif length(Neigh_ind) < 1
                    disp('H atom under coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'H_'};
                    if Heal_H == 1
                        atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    elseif length(Neigh_ind) < 1
                        Neighbours
                        atom(i)=[];
                        disp('removed atom...')
                        i
                        %pause
                    end
                    
                end
                
            end
            
            % If the element is O
            if strncmp(atom(i).type,{'O'},1)
                All_Neighbours=[All_Neighbours;Neighbours];
                if strcmp(Neighbours,'AlAlH')
                    atom(i).fftype={'Oh'};
                elseif strcmp(Neighbours,'AlAlSi')
                    atom(i).fftype={'Op'};
                elseif strncmp(Neighbours,'HHMgo',4)
                    atom(i).fftype={'Oalhh'};
                elseif strncmp(Neighbours,'HHMg',4)
                    atom(i).fftype={'Oalhh'};
                elseif strcmp(Neighbours,'AlHH')
                    atom(i).fftype={'Oalhh'};
                elseif strcmp(Neighbours,'AlHSi')
                    atom(i).fftype={'Oahs'}; % Al-OH-Si for acidic edge
                elseif strcmp(Neighbours,'AlH')
                    atom(i).fftype={'Oalh'}; % Al-O-H or Al-O-Si
                elseif  strncmp(Neighbours,'AlHMgo',5)
                    atom(i).fftype={'Ohmg'};
                elseif  strncmp(Neighbours,'AlHMg',5)
                    atom(i).fftype={'Ohmg'};
                elseif strcmp(Neighbours,'AlMgoSi')
                    atom(i).fftype={'Omg'};
                elseif strcmp(Neighbours,'AlMgSi')
                    atom(i).fftype={'Omg'};
                elseif strcmp(Neighbours,'HSi')
                    atom(i).fftype={'Osih'};
                elseif strcmp(Neighbours,'AlOmg')
                    atom(i).fftype={'Odsub'};
                elseif strcmp(Neighbours,'AlAltH')
                    atom(i).fftype={'Oalt'};
                elseif strcmp(Neighbours,'AlAlAlt')
                    atom(i).fftype={'Oalt'};
                elseif strcmp(Neighbours,'AlAlAl')
                    atom(i).fftype={'Oalt'};
                elseif strcmp(Neighbours,'AltSi')
                    atom(i).fftype={'Oalt'};
                elseif strcmp(Neighbours,'AlSi') % Al-O-H or Al-O-Si
                    if sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) < 1
                        atom(i).fftype={'Oalsi'};
                    elseif sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) > 0
                        atom(i).fftype={'Oalt'};
                    end
                    %                 elseif strcmp(Neighbours,'AlAl')
                    %                     atom(i).fftype={'O'};
                elseif strcmp(Neighbours,'FeFe') || strcmp(Neighbours,'FeoFeo')
                    atom(i).fftype={'Odsub'};
                elseif strcmp(Neighbours,'MgoMgoSi') || strcmp(Neighbours,'FeFeSi')
                    atom(i).fftype={'Odsub'};
                elseif strcmp(Neighbours,'MgMgSi') || strcmp(Neighbours,'FeoFeoSi')
                    atom(i).fftype={'Odsub'};
                elseif strcmp(Neighbours,'SiSi')
                    atom(i).fftype={'Ob'}; % Old version atom(i).fftype={'O'};
                elseif strcmp(Neighbours,'MnMnMn')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn4Mn4Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn3Mn4Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn2Mn3Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn2Mn4Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn4Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn3Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn2Mn3')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'Mn2Mn4')
                    atom(i).fftype={'Omn'};
                elseif strcmp(Neighbours,'AlAlSiSi')
                    atom(i).fftype={'Oz'};
                elseif strcmp(Neighbours,'HH')
                    atom(i).fftype={'Ow'};
                    disp('Water?')
                elseif length(Neigh_ind) > 2
                    disp('O atom over coordinated')
                    i
                    Neighbours
                    atom(i).fftype={'O^'};
                elseif length(Neigh_ind) == 1 || strcmp(Neighbours,'AlAl') || strcmp(Neighbours,'AlMg') || strcmp(Neighbours,'Si')
                    if strcmp(Neighbours,'Si')
                        atom(i).fftype={'Osi'};
                    elseif strncmp(Neighbours,'AlAl',2)
                        atom(i).fftype={'Oal'};
                    else
                        disp('O atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'O_'};
                    end
                    if Heal_O == 1
                        try
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                        catch
                            atom(size(atom,2)+1)=atom(i);
                        end
                        atom(end).type={'H'};
                        atom(end).fftype={'H'};
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                elseif length(Neigh_ind) == 0
                    Neighbours
                    atom(i)
%                     atom(i)=[];
                     disp('remove O atom...?')
                     
                     rm_ind=[rm_ind i];
                     
                end
                if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 %&& nH_extra > nH_added;
                    disp('Adding acidic H to Oalh-->Oalhh')
                    atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                    atom(end).index=size(atom,2);
                    Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                    NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -1/2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)));
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
                if strcmp(atom(i).fftype,'Oalsi') && Add_extra_H == 1 %&& nH_extra > nH_added;
                    disp('Adding acidic H to Oalsi-->Oahs')
                    atom(i).fftype={'Oahs'};
                    atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                    atom(end).index=size(atom,2);
                    Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                    NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)));
                    [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                end
            end
        end
        %         end
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

try
    % If the usual atom types
    atom = charge_atom(atom,Box_dim,ffname,watermodel);
    % If new atom types
    atom = charge_atom(atom,Box_dim,ffname,watermodel,'more');
catch
    disp('Could not set the charge...')
end
assignin('caller','newatom',atom);
assignin('caller','remove_ind',rm_ind);

for i=1:length(unique([atom.type]))
    try
        new_Atom_label=sort(unique([atom.type]));
        ind=ismember([atom.type],new_Atom_label(i));
        assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
    catch
        disp('Could not finalize:')
        i
        Atom_label
    end
end
ffname
watermodel

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),'SPC/E');
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));
% atom = charge_clayff_atom(atom,Box_dim)
% Total_charge

