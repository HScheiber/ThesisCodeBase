function Output_Structure = Shorten_Structure_Name(Input_Structure)
    Output_Structure = strrep(Input_Structure,'Rocksalt','Rock.');
    Output_Structure = strrep(Output_Structure,'Wurtzite','Wurtz.');
    Output_Structure = strrep(Output_Structure,'BetaBeO','$\beta$-BeO');
    Output_Structure = strrep(Output_Structure,'FiveFive','5-5');
    Output_Structure = strrep(Output_Structure,'Sphalerite','Sphal.');
end