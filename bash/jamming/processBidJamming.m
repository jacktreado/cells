function processBidJamming(simdir,savestr)
%% FUNCTION to loop over completed bid cell jamming simulations
%   -- simdir: parent directory with simulation directories as children
%   -- savestr: mat file string to store save data

% simdir string
if strcmp(simdir(end),'/')
    simdir(end) = [];
end

% print to console
fprintf('\t ** In processBidJamming, processing data from %s\n',simdir);

%% Get list of files

% get directory list
flist = dir([simdir '/bidcells*']);
NF = length(flist);
if NF == 0
    error('No directories found in processBidJamming with this input, ending.');
end


%% Loop over files, read data and extract

% data to save
phiJ = zeros(NF,1);
NvvJ = zeros(NF,1);
NccJ = zeros(NF,1);
calAJ = cell(NF,1);
calA0 = cell(NF,1);
vforceJ = cell(NF,2);
dskip = false(NF,1);

% loop over directories
for dd = 1:NF
    % directory info
    fname = flist(dd).name;
    jamfstr = [flist(dd).folder '/' fname];
    
    % check if jamfstr has any size
    finfo = dir(jamfstr);
    fsize = finfo.bytes;
    if fsize == 0
        fprintf('\t ** jam file %s empty, skipping...\n',fname);
        dskip(dd) = true;
        continue;
    end
    
    % load in data from jam string
    cellJamData = readCellJamData(jamfstr);
    
    % parse data for saving
    phiJ(dd)        = cellJamData.phi(1);
    NvvJ(dd)        = cellJamData.Nvv(1);
    NccJ(dd)        = cellJamData.Ncc(1);
    calA0{dd}       = cellJamData.calA0(1,:);
    calAJ{dd}       = cellJamData.calA(1,:);
    vforceJ{dd,1}   = cellJamData.xfrc(1,:);
    vforceJ{dd,2}   = cellJamData.yfrc(1,:);
end

% delete extra entries
phiJ(dskip) = [];
NvvJ(dskip) = [];
NccJ(dskip) = [];
calAJ(dskip) = [];
vforceJ(dskip,:) = [];

% save config-wide data
NV = cellJamData.nv(1,:);
NCELLS = cellJamData.NCELLS;
NCONFIGS = sum(~dskip);

%% save to savestr

save(savestr,'NCELLS','NCONFIGS','calA0','NV',...
    'phiJ','NvvJ','NccJ','calAJ','vforceJ','cellJamData');


end