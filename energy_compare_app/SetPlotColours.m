function plotObjs_sorted = SetPlotColours(plotObjs,Ctype,Color_Scheme)

% Remove empty
plotObjs = plotObjs(~cellfun(@isempty,plotObjs));

N = length(plotObjs);

% Pre allocate arrays
UserDataArray = cell(N,4);
SortedArray_OneLine = cell(N,1);

% Get user data
for i = 1:N
    UserDataArray(i,:) = plotObjs{i}.UserData;
end

% Sort by Basis set > Salt > Theory
[SortedArray,index] = sortrows(UserDataArray,[4 1 2 3]);
plotObjs_sorted = plotObjs(index);

% Ignore structure for deciding colours
SortedArray_NoStruc = SortedArray(:,[1:2 4]);

for i = 1:N
    SortedArray_OneLine{i} = [SortedArray_NoStruc{i,:}];
end

% Get number of colours
NCol = length(unique(SortedArray_OneLine));
if NCol < 3
    NCol = 3;
end

Colours = cbrewer(Ctype,Color_Scheme,NCol,'PCHIP');

% Assign colours and line/marker types
Previous = '';
q = 0;
for i = 1:N
    Current = SortedArray_OneLine{i};
    Structure = SortedArray{i,3};
    if ~strcmp(Previous,Current)
        q = q+1;
    end
    
    if strcmp(plotObjs_sorted{i}.Type,'line')
        plotObjs_sorted{i}.Color = Colours(q,:);
        
        if strcmp(Structure,'Rocksalt')
             plotObjs_sorted{i}.LineStyle = '-';
             plotObjs_sorted{i}.Marker = 'none';
        elseif strcmp(Structure,'Wurtzite')
             plotObjs_sorted{i}.LineStyle = '--';
             plotObjs_sorted{i}.Marker = 'none';
        elseif strcmp(Structure,'Pair')
             plotObjs_sorted{i}.LineStyle = ':';
             plotObjs_sorted{i}.Marker = 'none';
        elseif strcmp(Structure,'Sphalerite')
             plotObjs_sorted{i}.LineStyle = '-.';
             plotObjs_sorted{i}.Marker = '+';
        elseif strcmp(Structure,'CsCl')
             plotObjs_sorted{i}.LineStyle = '-';
             plotObjs_sorted{i}.Marker = 'x';
        elseif strcmp(Structure,'NiAs')
             plotObjs_sorted{i}.LineStyle = '--';
             plotObjs_sorted{i}.Marker = 's';
        elseif strcmp(Structure,'BetaBeO')
             plotObjs_sorted{i}.LineStyle = '-.';
             plotObjs_sorted{i}.Marker = '^';
        elseif strcmp(Structure,'FiveFive')
             plotObjs_sorted{i}.LineStyle = '-';
             plotObjs_sorted{i}.Marker = 'p';
        end
        
    elseif strcmp(plotObjs_sorted{i}.Type,'scatter')
        plotObjs_sorted{i}.CData = Colours(q,:);
        
        if strcmp(Structure,'Rocksalt')
            plotObjs_sorted{i}.Marker = 'o';
        elseif strcmp(Structure,'Wurtzite')
            plotObjs_sorted{i}.Marker = '*';
        elseif strcmp(Structure,'Pair')
            plotObjs_sorted{i}.Marker = '.';
        elseif strcmp(Structure,'Sphalerite')
            plotObjs_sorted{i}.Marker = '+';
        elseif strcmp(Structure,'CsCl')
            plotObjs_sorted{i}.Marker = 'x';
        elseif strcmp(Structure,'NiAs')
            plotObjs_sorted{i}.Marker = 's';
        elseif strcmp(Structure,'BetaBeO')
            plotObjs_sorted{i}.Marker = '^';
        elseif strcmp(Structure,'FiveFive')
            plotObjs_sorted{i}.Marker = 'p';
        end
        
    end
    Previous = Current;
end
end