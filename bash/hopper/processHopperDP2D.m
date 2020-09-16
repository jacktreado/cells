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
calA = cell(NF,1);
calA0 = cell(NF,1);

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
    NT              = dpHopperPosData.NFRAMES;
    NFRAMES(ff)     = NT;
    
    % -- get hopper geometry
    w               = dpHopperPosData.w;
    w0              = dpHopperPosData.w0;
    L               = dpHopperPosData.L;
    th              = dpHopperPosData.th;
    
    hopperInfo(ff,1) = w;
    hopperInfo(ff,2) = w0;
    hopperInfo(ff,3) = L;
    hopperInfo(ff,4) = th;
    
    % -- get particle info
    fprintf('\t ** -- Computing particle information\n');
    
    nv          = dpHopperPosData.nv;
    nv          = nv(1,:);
    xpos        = dpHopperPosData.xpos;
    ypos        = dpHopperPosData.ypos;
    a0          = dpHopperPosData.a0;
    l0          = dpHopperPosData.l0;
    
    
    % save shape data and com position data
    calA0tmp    = zeros(NT,N);
    calAtmp     = zeros(NT,N);
    cxpos       = zeros(NT,N);
    for tt = 1:NT
        for nn = 1:N
            % vertex positions
            x = xpos{tt,nn};
            y = ypos{tt,nn};
            
            % shape info
            nvtmp = nv(nn);
            a0tmp = a0(tt,nn);
            l0tmp = l0(tt,nn);
            
            % area
            atmp = polyarea(x,y);
            
            % perimeter
            lvx = x([2:end 1]) - x;
            lvy = y([2:end 1]) - y;
            l = sqrt(lvx.^2 + lvy.^2);
            ptmp = sum(l);
            
            % shape parameters
            calAtmp(tt,nn) = ptmp^2/(4.0*pi*atmp);
            calA0tmp(tt,nn) = (nvtmp*l0tmp)^2/(4.0*pi*a0tmp);
            
            % cxpos
            cxpos(tt,nn) = mean(x);
        end
    end
    calA{ff} = calAtmp;
    calA0{ff} = calA0tmp;
    
    % Loop over time, find outflow events
    fprintf('\t ** -- Computing outflow statistics\n');
    dnout = zeros(NT-1,1);
    for tt = 2:NT
        flowedOut = (cxpos(tt,:) > L & cxpos(tt-1,:) < L);
        dnout(tt-1) = sum(flowedOut);
    end
    nout{ff} = cumsum(dnout);
    
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
wallFrc(rmvf,:)     = [];
nout(rmvf,:)        = [];

%% Save to matfile
save(savestr,...
    'wallFrc','nout','calA0','calA','hopperInfo','NFRAMES');

end






