function processHopperDP2D(N,simstrpattern,savestr)
%% Load in hopper data, extract data, save to save string
% output statistics:
%   -- net forces on walls
%   -- outflow statistics (nout)
%   -- shape data
%   -- speed data

% get files
flist   = dir([simstrpattern '*']);
NF      = length(flist);
if NF == 0
    fprintf('\t ** ERROR: no files found for simstrpattern = %s, ending.\n',simstrpattern);
    error('No files found');
end

%% Loop over files, save data from each sim

% data arrays

% remover
rmvf = false(NF,1);

% simulation info
NFRAMES = zeros(NF,1);
hopperInfo = zeros(NF,4);

% particle statistics
calA = cell(NF,N);
calA0 = cell(NF,N);
speed = cell(NF,N);

% wall forces
wallFrc = cell(NF,4);

% outflow statistics
nout = cell(NF,1);

for ff = 1:NF
    % get file information
    fname = flist(ff).name;
    ffldr = flist(ff).folder;
    fstr = [ffldr '/' fname];
    
    % check file size
    finfo = dir(fstr);
    if finfo.bytes == 0
        fprintf('\t** File %s empty skipping.\n',fname);
        rmvf(ff) = true;
        continue;
    else
        fprintf('\t** Reading in %s\n',fname);
    end
    
    % read in data
    dpHopperPosData = readDPHopperData(fstr);
    
    % -- get number of frames in simulation
    NF              = dpHopperPosData.NF;
    NFRAMES(ff)     = NF;
    
    % -- get hopper geometry
    w               = dpHopperPosData.w;
    w0              = dpHopperPosData.w0;
    L               = dpHopperPosData.L;
    th              = dpHopperPosData.th;
    
    hopperInfo(ff,1) = w;
    hopperInfo(ff,2) = w0;
    hopperInfo(ff,3) = L;
    hopperInfo(ff,4) = th;
    
    % -- get outflow statistics
    fprintf('\t ** -- Computing outflow statistics\n');
    cxpos = dpHopperPosData.cxpos;
    
    % Loop over time, find outflow events
    dnout = zeros(NF-1,1);
    for tt = 2:NF
        flowedOut = (cxpos(tt,:) > L & cxpos(tt-1,:) < L);
        dnout(tt-1) = sum(flowedOut);
    end
    nout{ff} = cumsum(dnout);
    
    % -- get particle info
    fprintf('\t ** -- Computing particle information\n');
    nv          = dpHopperPosData.nv;
    nv          = nv(1,:);
    calAtmp     = dpHopperPosData.calA;
    calA0tmp    = dpHopperPosData.calA0;
    xvel        = dpHopperPosData.xvel;
    yvel        = dpHopperPosData.yvel;
    for nn = 1:N
        % get shapes
        calA{ff,nn}     = calAtmp(:,nn);
        calA0{ff,nn}    = calA0tmp(:,nn);
        
        % get speeds
        speednn = zeros(NT,1);
        for tt = 1:NT
            vcomx           = sum(xvel{tt,nn})/nv(nn);
            vcomy           = sum(yvel{tt,nn})/nv(nn);
            speednn(tt)     = sqrt(vcomx^2 + vcomy^2);
        end
        
        % save speed
        speed{ff,nn}    = speednn;
    end
    
    % -- get wall forces
    fprintf('\t ** -- Saving wall forces\n');
    wftmp           = dpHopperPosData.wallfrc;
    wallFrc{ff,1}   = wftmp(:,1);
    wallFrc{ff,2}   = wftmp(:,2);
    wallFrc{ff,3}   = wftmp(:,3);
    wallFrc{ff,4}   = wftmp(:,4);
    
    
    fprintf('\t ** -- Finished parsing data for %s\n\n',fname);
end

% remove extra slots
NFRAMES(rmvf)       = [];
hopperInfo(rmvf,:)  = [];
calA(rmvf,:)        = [];
calA0(rmvf,:)       = [];
speed(rmvf,:)       = [];
wallFrc(rmvf,:)     = [];
nout(rmvf,:)        = [];

%% Save to matfile
save(savestr,...
    'wallFrc','speed','nout','calA0','calA','hopperInfo','NFRAMES');
