function cellTrajectoryData = readQuenchedPackings(fstr)
%% FUNCTION to read in cell trajectoryu data given file string

% print info to console
finfo = dir(fstr);
fprintf('Reading in %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% open file stream
fid = fopen(fstr);

% read in sim details from first frame
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);   NCELLS = NCELLS{1};
phi0        = textscan(fid,'PACKF %f',1);                   phi0 = phi0{1};
fline       = fgetl(fid);
Ltmp        = textscan(fid,'BOXSZ %f %f',1);
fline       = fgetl(fid);

% cells to save 
NFRAMES = 5e4;
nv      = zeros(NFRAMES,NCELLS);
xpos    = cell(NFRAMES,NCELLS);
ypos    = cell(NFRAMES,NCELLS);
zc      = zeros(NFRAMES,NCELLS);
zv      = zeros(NFRAMES,NCELLS);
l0      = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
L       = zeros(NFRAMES,2);

% number of frames found
nf = 1;

% loop over frames, read in data
while ~feof(fid)
    % get box length
    L(nf,1) = Ltmp{1};
    L(nf,2) = Ltmp{2};
    
    % get info about deformable particle
    for nn = 1:NCELLS
        % get cell pos and asphericity
        cInfoTmp = textscan(fid,'CINFO %f %f %f %f %f',1);   
        fline = fgetl(fid);     % goes to next line in file

        NVTMP = cInfoTmp{1};
        nv(nf,nn) = NVTMP;
        zc(nf,nn) = cInfoTmp{2};
        zv(nf,nn) = cInfoTmp{3};
        a0(nf,nn) = cInfoTmp{4};
        l0(nf,nn) = cInfoTmp{5};
        
        % get vertex positions
        vPosTmp = textscan(fid,'VINFO %*f %*f %f %f',NVTMP); 
        fline = fgetl(fid);     % goes to next line in file

        % parse data
        xposTmp = vPosTmp{1};
        yposTmp = vPosTmp{2};

        % save in cell
        xpos{nf,nn} = xposTmp;
        ypos{nf,nn} = yposTmp;
    end
    
    % increment frame count
    nf = nf + 1;
    
    % read in/discard trash information
    % NOTE: IF HOPPER GEOMETRY/NUMBER OF PARTICLES/ANYTHING ELSE CHANGES
    % DURING SIMULATION, THIS IS WHERE YOU WOULD CHECK AND MAKE ADJUSTMENTS
    
    % get ENDFR string, check that read is correct
    fline = fgetl(fid);
    endFrameStr = sscanf(fline,'%s');
    
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
        
        % read in sim details from first frame
        NCELLStmp       = textscan(fid,'NUMCL %f');
        
        % read in packing fraction
        phi0            = textscan(fid,'PACKF %f',1);                   
        phi0            = phi0{1};
        fline = fgetl(fid);
        
        % update box size
        Ltmp            = textscan(fid,'BOXSZ %f %f',1);
        L(nf,1)         = Ltmp{1};
        L(nf,2)         = Ltmp{2};
        fline = fgetl(fid);
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    xpos(nf:end,:) = [];
    ypos(nf:end,:) = [];
    nv(nf:end,:) = [];
    zc(nf:end,:) = [];
    zv(nf:end,:) = [];
    l0(nf:end,:) = [];
    a0(nf:end,:) = [];
    L(nf:end,:) = [];
end

% close position file
fclose(fid);

% store cell pos data into struct
cellTrajectoryData              = struct('NFRAMES',NFRAMES,'NCELLS',NCELLS);
cellTrajectoryData.phi0         = phi0;
cellTrajectoryData.xpos         = xpos;
cellTrajectoryData.ypos         = ypos;
cellTrajectoryData.nv           = nv;
cellTrajectoryData.zc           = zc;
cellTrajectoryData.zv           = zv;
cellTrajectoryData.l0           = l0;
cellTrajectoryData.a0           = a0;
cellTrajectoryData.L            = L;

end