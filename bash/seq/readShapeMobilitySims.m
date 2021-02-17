function readShapeMobilitySims(da,savestr)
%% FUNCTION to read through all shape mobility sims
% extract parameters, save long-time shape fluctuations

% get list of directories
dirlist = dir(datalocinfo);
ND = length(dirlist);

% parameter lists
Kb = zeros(ND,1);
v0 = zeros(ND,1);


% loop over directories

end