function map_out = center_colormap(map,data,indexValue)

% Calculate where proportionally indexValue lies between minimum and
% maximum values
L = length(data);
maxdat = max(data);
mindat = min(data);
index = L*abs(indexValue-mindat)/(maxdat-mindat);

% Get map top, middle, and bottom colors
N = size(map,1);
Mid_point = ceil(N/2);
topColors = map(Mid_point+1:end,:);
indexColor = map(Mid_point,:);
bottomColors = map(1:Mid_point-1,:);

% Loop through bottom colors
cmap = [];
for idx = 1:size(bottomColors,1)-1
    customCMap = [linspace(bottomColors(idx,1),bottomColors(idx+1,1),5*index)',...
            linspace(bottomColors(idx,2),bottomColors(idx+1,2),5*index)',...
            linspace(bottomColors(idx,3),bottomColors(idx+1,3),5*index)'];
    cmap = [cmap; customCMap];
end

% Do the middle ones
customCMap = [linspace(bottomColors(end,1),indexColor(1),5*index)',...
        linspace(bottomColors(end,2),indexColor(2),5*index)',...
        linspace(bottomColors(end,3),indexColor(3),5*index)'];
cmap = [cmap; customCMap];

customCMap = [linspace(indexColor(1),topColors(1,1),5*(L-index))',...
        linspace(indexColor(2),topColors(1,2),5*(L-index))',...
        linspace(indexColor(3),topColors(1,3),5*(L-index))'];
cmap = [cmap; customCMap];


% Do the top colors
for idx = 1:size(topColors,1)-1
    customCMap = [linspace(topColors(idx,1),topColors(idx+1,1),5*(L-index))',...
            linspace(topColors(idx,2),topColors(idx+1,2),5*(L-index))',...
            linspace(topColors(idx,3),topColors(idx+1,3),5*(L-index))'];
    cmap = [cmap; customCMap];
end

map_out = cmap;  % Combine colormaps

end