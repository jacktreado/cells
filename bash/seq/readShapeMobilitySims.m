function readShapeMobilitySims(da,savestr)
%% FUNCTION to read through all shape mobility sims
% extract parameters, save long-time shape fluctuations

% get list of directories
dirlist = dir(datalocinfo);
ND = length(dirlist);

% parameter lists
Kb = zeros(ND,1);
v0 = zeros(ND,1);
Dr = zeros(ND,1);
NT = zeros(ND,1);

% data
meanCalA = zeros(ND,1);


% loop over directories, extract parameters
for dd = 1:ND

end