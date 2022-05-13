%% analyze_atom.m
% * This function fetches the ionic radius, originally taken from below
% * Ref. 1	Revised effective ionic radii and systematic studies of
% interatomic distances in halides and chalcogenides. R. D. Shannon Acta
% Cryst. (1976) A32, 751-767.
% * Ref. 2	Electronic Table of Shannon Ionic Radii, J. David Van Horn,
% 2001, downloaded MO/DA/YEAR.
% *
% * This function also calculates the bond valence values according to
% http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% compiled by I. David Brown, McMaster University, Ontario, Canada
% * Data set bvparm2016.cif: 2016 version, (posted 2016-11-03)
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # analyze = analyze_atom(atom,Box_dim)
% # analyze = analyze_atom(atom,Box_dim,2.5)

function properties = analyze_atom(atom,Box_dim,varargin)

% FIX multiple OxStates and other from the Shannon file

load('Revised_Shannon_radii.mat');

if nargin > 2
    rmax=varargin{1};
else
    rmax=2.25;
end

element=element_atom(atom); % Optionally, set the names ow water O and water H
%element=bond_angle_atom(element,Box_dim,1.25,rmax,'more');
disp('Trying to find bonded atoms')
element=bond_atom(element,Box_dim,rmax);
disp('Running Bond Valence analysis...')
element=bond_valence_atom(element,Box_dim,1.25,rmax);

% other_ind=[];
% for i=1:size(element,2)
%     try
%         ValenceGuess=round([element(i).valence]);
%         element(i)=bond_valence_atom(element(i),Box_dim,1.25,rmax,ValenceGuess);
%     catch
%         other_ind=[other_ind i];
%     end
%     if mod(i,100)==1
%         i-1
%     end
% end
% if numel(other_ind)>0
%     element(other_ind)=bond_valence_atom(element(other_ind),Box_dim,1.25,rmax);
% end

element=mass_atom(element);
assignin('caller','element',element);
for i=1:size(element,2)
    Ion_ind=find(strcmp([element(i).type],Ion));
    if numel(Ion_ind)==0
        Ion_ind=find(strncmpi([element(i).type],Ion,1));
    end
    CN_ind=find(numel(element(i).neigh.index)==CN);
    Ox_ind=find(round(element(i).valence)==OxState);
    CN_ind=intersect(Ion_ind,CN_ind);
    Ox_ind=intersect(Ion_ind,Ox_ind);
    ind=intersect(Ox_ind,CN_ind);
    if numel(ind)==0
        ind=Ion_ind(1);
    elseif numel(ind)>1
        ind=ind(1);
    end
    properties(i).index=element(i).index;
    properties(i).type=Ion(ind);
    properties(i).fftype=atom(i).type;
    properties(i).neigh=element(i).neigh;
    properties(i).bond=element(i).bond;
    %     properties(i).angle=element(i).angle; % Because bond_atom() and not bond_angle_atom()
    properties(i).bv=element(i).bv;
    properties(i).valence=element(i).valence;
    properties(i).ave_dist=mean(properties(i).bond.dist);
    properties(i).std_dist=std(properties(i).bond.dist);
    properties(i).rdiffvalence=element(i).Rdiff;
    properties(i).cn_bv=size(properties(i).bv,2);
    properties(i).ShannonParam={'-->'};
    properties(i).atnum=element(i).atnum;
    properties(i).mass=element(i).mass;
    properties(i).oxstate=OxState(ind);
    properties(i).cn_guess=CN(ind);
    properties(i).ionicradii=IonicRadii(ind);
    if mod(i,100)==1
        i-1
    end
end

assignin('caller','pre_properties',properties);

for i=1:size(properties,2)
    Ion_ind=find(strcmp([element(i).type],Ion));
    if numel(Ion_ind)==0
        Ion_ind=find(strncmpi([element(i).type],Ion,1));
    end
    CN_ind=find(numel(element(i).neigh.index)==CN);
    ind=intersect(Ion_ind,CN_ind);
    
    if numel(ind)>0
        current_radii=[properties(i).ionicradii];
        neigh_radii=[properties([properties(i).neigh.index]).ionicradii];
        sum_radii=repmat(current_radii',numel(neigh_radii),1)+repmat(neigh_radii,numel(current_radii),1)';
        [minvalue,preox_ind]=min(abs(mean(sum_radii-mean([properties(i).neigh.dist]-properties(i).rdiffvalence))));
        preox_ind;
        properties(i).type=Ion(ind(preox_ind));
%         properties(i).oxstate=OxState(ind(preox_ind));
        properties(i).ip=ZoverIR(ind(preox_ind));
        properties(i).cn_guess=CN(ind(preox_ind));
        properties(i).crysradii=CrysRadii(ind(preox_ind));
        properties(i).ionicradii=IonicRadii(ind(preox_ind));
        properties(i).vdwradii=radius_vdw([element(i).type]);
        properties(i).ip=ZoverIR(preox_ind);
        properties(i).elecconf=strtrim(ElecConf(ind(preox_ind)));
        properties(i).spinstate=strtrim(SpinState(ind(preox_ind)));
    end
    if mod(i,100)==1
        i-1
    end
end

diff_valence=(abs([properties.oxstate])-[properties.valence])';
ind=find(abs(diff_valence)>0.5);
if numel(ind)>0
    disp('Possible problems with index due to multiple oxidation states,')
    disp('or atoms that are undersaturated!')
    unique([atom(ind).type])
    ind
    assignin('caller','heal_ind',ind');
end
assignin('caller','heal_ind',ind');

disp('Global instability index is:')
GII=(sum((abs([properties.oxstate])-[properties.valence]).^2)/size(properties,2))^0.5
if GII>0.2
    disp('GII > 0.2 --> Structure likely not super stable...');
end

Atom_labels=unique([atom.fftype]);
for i=1:length(Atom_labels)
    ind=find(strcmp([atom.fftype],Atom_labels(i)));
    BondSummary(i).type=Atom_labels(i);
    BondSummary(i).GII=(sum((abs([properties(ind).oxstate])-[properties(ind).valence]).^2)/numel(ind))^0.5;
    BondSummary(i).d_strain=sum([properties(ind).valence]-abs([properties(ind).oxstate]))/numel(ind);
    BondSummary(i).ValenceAve=mean([properties(ind).valence]);
    BondSummary(i).ValenceStd=std([properties(ind).valence]);
    BondSummary(i).dist=mean([properties(ind).ave_dist]);
end

load('bond_valence_values.mat');
for i=1:size(Bond_index,1)
    Bond_index(i,4)=properties(Bond_index(i,1)).ionicradii+properties(Bond_index(i,2)).ionicradii;
    [mean_bv,std_bv,bv,bvalue]=bond_valence_data(properties(Bond_index(i,1)).type,properties(Bond_index(i,2)).type,Bond_index(i,3),Ion_1,Ion_2,R0,b,Valence_1,Valence_2);
    Bond_index(i,5)=bv;
end

format short
diff_bond=Bond_index(:,4)-Bond_index(:,3);
ind=find(abs(diff_bond)>0.5);
if numel(ind)>0
    disp('Possible problems with bond between:')
    for i=1:numel(ind)
        [Bond_index(ind(i),1) Bond_index(ind(i),2)]
        properties(Bond_index(ind(i),1)).type
        properties(Bond_index(ind(i),2)).type
        Bond_index(ind(i),3:end)
    end
end

assignin('caller','GII',GII);
assignin('caller','BondSummary',BondSummary);
assignin('caller','diff_bond',diff_bond);
assignin('caller','diff_valence',diff_valence);
assignin('caller','prop_atom',element);
try
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','diff_bond_bv',[properties.rdiffvalence]');
catch
    disp('Could not assignin Bond_index, Angle_index, distmatrix')
end
