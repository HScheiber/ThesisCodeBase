% Function that force closes all figures, even those with modified close
% request functions
figHandles = findobj('Type', 'figure');
for i=1:length(figHandles)
    Currentfig = figHandles(i);
    if isvalid(Currentfig)
        Currentfig.CloseRequestFcn = 'closereq';
    end
end

close(figHandles);