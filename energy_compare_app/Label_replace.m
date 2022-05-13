function labelText = Label_replace(Label)
    labelText = strrep(Label,'Rocksalt','_R');
    labelText = strrep(labelText,'Wurtzite','_W');
    labelText = strrep(labelText,'Sphalerite','_S');
    labelText = strrep(labelText,'CsCl','_C');
    labelText = strrep(labelText,'NiAs','_N');
    labelText = strrep(labelText,'BetaBeO','_B');
    labelText = strrep(labelText,'FiveFive','_F');
    labelText = strrep(labelText,'Pair','_P');
end