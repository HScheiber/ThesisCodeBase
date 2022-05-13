function parallelsave(fname____,vnames___,varargin)
  numvars___=numel(varargin);
  for ii___=1:numvars___
       eval([vnames___{ii___},'=varargin{ii___};']);  
  end
  save('-mat',fname____,vnames___{1});
  for ii___ = 2:numvars___    
    save('-mat',fname____,vnames___{ii___},'-append');
  end
end