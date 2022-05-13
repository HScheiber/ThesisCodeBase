function [Bar_Y_Trimmed,Bar_SD_Trimmed,BarCData,Bar_Label_Trimmed] = SetPlotColoursBar(Bar_Y,Bar_Label,Bar_Info,Ctype,Color_Scheme)

% Remove empty entries
Empty_Index = ~arrayfun(@isnan,Bar_Y);
Bar_Y = Bar_Y(Empty_Index);
Bar_Label = Bar_Label(Empty_Index);

N = length(Bar_Y);

for i=length(Bar_Info):-1:1
    if ~Empty_Index(i)
        Bar_Info(i,:) = [];
    end
end

% Sort by Structure > Salt > Basis set > Theory
[Bar_Info_Sorted,index] = sortrows(Bar_Info,[3 1 4 2]);
Structures_Sorted = Bar_Info(index,3);
Bar_Y_Sorted = Bar_Y(index);
Bar_Label_Sorted = Bar_Label(index);

% Average experimental values
for i = 1:N
    Bar_Info_Sorted_OneLine{i} = [Bar_Info_Sorted{i,:}];
end
[~,ia,ic] = unique(Bar_Info_Sorted_OneLine,'stable');
N = length(ia);
Structures_Trimmed = Structures_Sorted(ia);
Bar_Label_Trimmed = Bar_Label_Sorted(ia);

Bar_Y_Trimmed = nan(N,1);
Bar_SD_Trimmed = nan(N,1);
for i = 1:N
    Bar_Y_Trimmed(i) = mean(Bar_Y_Sorted(ic == i));
    Bar_SD_Trimmed(i) = std(Bar_Y_Sorted(ic == i));
end

% Get number of colours
NCol = length(unique(Structures_Trimmed));
if NCol < 3
    NCol = 3;
end

Colours = cbrewer(Ctype,Color_Scheme,NCol,'PCHIP');

% Assign colours
Previous = '';
q = 0;
BarCData = nan(N,3);
for i = 1:N
    Current = Structures_Trimmed{i};
    if ~strcmp(Previous,Current)
        q = q+1;
    end
    
    BarCData(i,:) = Colours(q,:);
    Previous = Current;
end

end