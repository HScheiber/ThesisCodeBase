function num = atomlabel2num(label)
if strcmp(label,'LI')
    num = '3';
elseif strcmpi(label,'NA')
    num = '11';
elseif strcmpi(label,'F')
    num = '9';
elseif strcmpi(label,'Cl')
    num = '17';
elseif strcmpi(label,'Br')
    num = '35';
elseif strcmpi(label,'I')
    num = '253';
end