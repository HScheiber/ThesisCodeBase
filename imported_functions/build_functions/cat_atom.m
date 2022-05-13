%% cat_atom.m
% * This is a special script (and not a function) that imports and appends atom structs into a
% * .gro trajectory file. It is useful for generating trajectories where
% * water molecules are evaporated during the simulations, since VMD and
% * other trajectory viewers cannot handle non-constant number of particles
% * in a trajectory

atom0=import_atom(strcat('evap_0.gro'));
% traj=zeros(125,3*size(atom,2));
% frame=zeros(1,3*size(atom,2));
% All_Box_dim=zeros(125,9);
% frames1=[0:20:120];frames2=[121:240];
% frames=sort([frames1 frames2]);

traj=zeros(120,3*size(atom,2));
frame=zeros(1,3*size(atom,2));
All_Box_dim=zeros(120,9);
frames=[0:1:120];

for i=1:numel(frames)
   i
   j=frames(i)
    try
      atom=import_atom(strcat('evap_',num2str(j),'.gro'));
      tempnum=3*size(atom,2);
      Xdata=XYZ_data(:,1);
      Ydata=XYZ_data(:,2);
      Zdata=XYZ_data(:,3);% -XYZ_data(1,3)+3.75;
%       Zdata(Zdata>Box_dim(3))=Zdata(Zdata>Box_dim(3))-Box_dim(3);
%       Zdata(Zdata<0)=Zdata(Zdata<0)+Box_dim(3);
      frame(1:3:tempnum)=Xdata; %0;
      frame(2:3:tempnum)=Ydata;
      frame(3:3:tempnum)=Zdata;
      traj(j+1,1:tempnum)=frame(1,1:tempnum);
      All_Box_dim(j+1,1:numel(Box_dim))=Box_dim;
   catch
       disp('No frame')
   end
end

traj=traj(frames+1,:);
Box_dim=All_Box_dim(frames+1,:);
if sum(Box_dim(:,9))==0
    Box_dim=Box_dim(:,1:3);
else
    Box_dim=Box_dim;
end

write_gro_traj(atom0,traj,Box_dim,'3D_all.gro')