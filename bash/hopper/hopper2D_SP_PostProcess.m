function hopper2D_SP_PostProcess(fileFolder,fileName,saveFile)
%% FUNCTION to parse multiple runs of 2D soft-particle hopper flows
% -- INPUT VARIABLES
%   ** fileFolder: directory of subdirectories with simulation data
%       -- NOTE: assume subdirs are of form [str]_seed*
% 
%   ** saveFile: .mat file, save position, vel, stress data
% 
% 

%% Get list of subdirectories

% check that fileFolder string does not end in an /
if strcmp(fileFolder(end),'/')
    fprintf('** Editing fileFolder string from %s ',fileFolder);
    fileFolder(end) = [];
    fprintf('to %s\n',fileFolder);
end

% list all subdirectories
subdirList = dir([fileFolder '/' fileName '_seed*']);

% number of subdirs
NSUBDIR = length(subdirList);
if NSUBDIR == 0
    error('No subdirectories found in this fileFolder, ending');
end

%% Loop over subdirectories, parse individual position files

% data to save
radiiList       = cell(NSUBDIR,1);  % particle radii
xposList        = cell(NSUBDIR,1);  % x coordinates
yposList        = cell(NSUBDIR,1);  % y coordinates
xvelList        = cell(NSUBDIR,1);  % x velocities
yvelList        = cell(NSUBDIR,1);  % y velocities
xfrcList        = cell(NSUBDIR,1);  % x forces
yfrcList        = cell(NSUBDIR,1);  % y forces
virialList      = cell(NSUBDIR,1);  % virial stresses
outflowList     = cell(NSUBDIR,1);  % outflow over time
hopperPhiList   = cell(NSUBDIR,1);  % hopper packing fraction over time

% save hopper info (w0, w, th)
hopperInfo = zeros(NSUBDIR,3);

% save number of particles
N           = zeros(NSUBDIR,1);

% save seeds
seeds       = zeros(NSUBDIR,1);

% loop over subdirectories, parse data in *.pos files,
% calc outflow and hopper packing fractions
for ss = 1:NSUBDIR
    
    % get file string
    fname = [subdirList(ss).name '.pos'];
    fstr = [subdirList(ss).folder '/' subdirList(ss).name '/' fname];
    
    % get seed
    seeds(ss) = sscanf(subdirList(ss).name,[fileName '_seed%d']);
    
    % read in data from pos file
    fid = fopen(fstr);

    % read in sim details from first frame
    NCELLS              = textscan(fid,'NUMCL %f',1,'HeaderLines',1);    
    NCELLS              = NCELLS{1};
    N(ss)               = NCELLS;
    hopperData          = textscan(fid,'HOPPR %f %f %f',1);
    hopperInfo(ss,1)    = hopperData{1};
    hopperInfo(ss,2)    = hopperData{2};
    hopperInfo(ss,3)    = hopperData{3};
    
    % store values
    w0  = hopperInfo(ss,1);
    w   = hopperInfo(ss,2);
    th  = hopperInfo(ss,3);
    
    % initialize arrays to store per frame data
    NFRAMES = 1e5;
    rad     = zeros(NFRAMES,NCELLS);
    xpos    = zeros(NFRAMES,NCELLS);
    ypos    = zeros(NFRAMES,NCELLS);
    xvel    = zeros(NFRAMES,NCELLS);
    yvel    = zeros(NFRAMES,NCELLS);
    xfrc    = zeros(NFRAMES,NCELLS);
    yfrc    = zeros(NFRAMES,NCELLS);
    virial  = zeros(NFRAMES,4);
    
    % number of frames found
    nf = 1;
    
    % loop over lines of file, parse frame by frame
    while ~feof(fid)
        % print to console
        fprintf('On frame %d of file %s\n',nf,fname);
        
        % get virial stress information
        vstress = textscan(fid,'VRIAL %f %f %f %f',1);
        virial(nf,1) = vstress{1};
        virial(nf,2) = vstress{2};
        virial(nf,3) = vstress{3};
        virial(nf,4) = vstress{4};

        % discard header line
        fline = fgetl(fid);
        fline = fgetl(fid);

        % get info about soft particle
        for ii = 1:NCELLS
            % read in data from file
            data        = textscan(fid,'SCELL %*d %f %f %f %f %f %f %f',1);

            % parse data
            rad(nf,ii)  = data{1};
            xpos(nf,ii) = data{2};
            ypos(nf,ii) = data{3};
            xvel(nf,ii) = data{4};
            yvel(nf,ii) = data{5};
            xfrc(nf,ii) = data{6};
            yfrc(nf,ii) = data{7};
        end

        % increment frame count
        nf = nf + 1;

        % read in/discard trash information
        % NOTE: IF HOPPER GEOMETRY/NUMBER OF PARTICLES/ANYTHING ELSE CHANGES
        % DURING SIMULATION, THIS IS WHERE YOU WOULD CHECK AND MAKE ADJUSTMENTS

        % get ENDFR string, check that read is correct
        fline = fgetl(fid);
        fline = fgetl(fid);
        endFrameStr = sscanf(fline,'%s %*s');

        if ~strcmp(endFrameStr,'ENDFR')
            error('Miscounted number of particles in .pos file, ending data read in');
        end

        % start new frame information
        if ~feof(fid)
            % get NEWFR line
            fline       = fgetl(fid);
            newFrameStr = sscanf(fline,'%s %*s');
            if ~strcmp(newFrameStr,'NEWFR')
                error('NEWFR not encountered when expected when heading to new frame...check line counting. ending.');
            end

            % read in sim details
            NCELLSdata  = textscan(fid,'NUMCL %f',1,'HeaderLines',1);    
            NCELLStmp   = NCELLSdata{1};
            hopperData  = textscan(fid,'HOPPR %f %f %f',1);
        end
    end

    % delete extra information
    if (nf < NFRAMES)
        NFRAMES = nf-1;
        rad(nf:end,:) = [];
        xpos(nf:end,:) = [];
        ypos(nf:end,:) = [];
        xvel(nf:end,:) = [];
        yvel(nf:end,:) = [];
        xfrc(nf:end,:) = [];
        yfrc(nf:end,:) = [];
        virial(nf:end,:) = [];
    end

    % calculate length
    L = 0.5*(w0 - w)/tan(th);

    % close position file
    fclose(fid);
    
    % store parsed variables to cells
    radiiList{ss}   = rad;
    xposList{ss}    = xpos;
    yposList{ss}    = ypos;
    xvelList{ss}    = xvel;
    yvelList{ss}    = yvel;
    xfrcList{ss}    = xfrc;
    yfrcList{ss}    = yfrc;
    virialList{ss}  = virial;    
    
% --- calculate outflow
    fprintf('\t -- saving outflow...\n');
    
    % number flowing out during a single frame
    instantOutflow = zeros(NFRAMES,1);

    % loop over frames, count number of x values that are greater than L
    for ff = 1:NFRAMES
        % all x values
        xtmp = xpos(ff,:);

        % add to outflow
        instantOutflow(ff) = sum(xtmp > L);
    end

    % get cumulative sum
    outflow = cumsum(instantOutflow);
    
    % save to list
    outflowList{ss} = outflow;
    
% --- calculate hopper packing fraction
    fprintf('\t -- saving hopper packing fraction...\n');

    % packing fraction over time
    hopperPhi = zeros(NFRAMES,1);
    
    % loop over frames, get phi approx based on reservoir or xmin
    for ff = 1:NFRAMES
        % all x values
        xtmp = xpos(ff,:);
        
        % get radii
        rtmp = rad(ff,:);

        % get min
        xmin = min(xtmp);
        
        % check which Lmin to use
        if xmin < -L
            Lmin = L;
        else
            if xmin < 0
                Lmin = abs(xmin);
            else
                Lmin = 0;
            end
        end
        
        % get hopper area
        aR = w0*Lmin;
        aH = L*0.5*(w0 + w);
        aT = aR + aH;
        
        % get approx area of disks
        hopperDisks     = xtmp > -L;
        hopperDiskRadii = rtmp(hopperDisks);
        diskArea        = pi*sum(hopperDiskRadii.^2);
        
        % save to time dependent packing fraction
        hopperPhi(ff)   = diskArea/aT;
    end
    
    % save to list
    hopperPhiList{ss} = hopperPhi;
end

% save to matfile
save(saveFile,'radiiList','xposList','yposList','xvelList','yvelList',...
    'xfrcList','yfrcList','virialList','outflowList','hopperPhiList',...
    'hopperInfo','N','seeds','subdirList');



end