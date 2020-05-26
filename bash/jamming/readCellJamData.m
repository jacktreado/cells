function cellJamData = readCellJamData(fstr)
%% FUNCTION to read in cell position data given file string

% print info to console
finfo = dir(fstr);
fprintf('Reading in %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% open file stream
fid = fopen(fstr);

% read in sim details from first frame
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);   NCELLS = NCELLS{1};
phi0        = textscan(fid,'PACKF %f',1);                   phi0 = phi0{1};
fline = fgetl(fid);
Ltmp        = textscan(fid,'BOXSZ %f %f',1);
fline = fgetl(fid);

% get contact info
cijtmp      = textscan(fid,'NCONTS %f %f',1);

% initialize arrays to store per frame data
NFRAMES = 5e4;
phi     = zeros(NFRAMES,1);
Ncc     = zeros(NFRAMES,1);
Nvv     = zeros(NFRAMES,1);
xpos    = cell(NFRAMES,NCELLS);
ypos    = cell(NFRAMES,NCELLS);
xvel    = cell(NFRAMES,NCELLS);
yvel    = cell(NFRAMES,NCELLS);
xfrc    = cell(NFRAMES,NCELLS);
yfrc    = cell(NFRAMES,NCELLS);
nv      = zeros(NFRAMES,NCELLS);
cxpos   = zeros(NFRAMES,NCELLS);
cypos   = zeros(NFRAMES,NCELLS);
l0      = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
del     = zeros(NFRAMES,NCELLS);
calA    = zeros(NFRAMES,NCELLS);
calA0   = zeros(NFRAMES,NCELLS);
L       = zeros(NFRAMES,2);

% number of frames found
nf = 1;

% loop over lines of file, parse frame by frame
while ~feof(fid)
    % get box length
    L(nf,1) = Ltmp{1};
    L(nf,2) = Ltmp{2};
    
    % store packing fraction and contacts 
    phi(nf) = phi0;
    Ncc(nf) = cijtmp{1};
    Nvv(nf) = cijtmp{2};
    
    % discard header line
    fline = fgetl(fid);
    
    % get info about deformable particle
    for nn = 1:NCELLS
        nvertTmp = textscan(fid,'NVERT %f',1); fline = fgetl(fid);     % goes to next line in file
        NVERT = nvertTmp{1};
        nv(nf,nn) = NVERT;
       
        % get cell pos and asphericity
        cInfoTmp = textscan(fid,'CELLP %*f %f %f %f %f %f %f %*f %f',1);   fline = fgetl(fid);     % goes to next line in file
        cxpos(nf,nn) = cInfoTmp{1};
        cypos(nf,nn) = cInfoTmp{2};
        l0(nf,nn) = cInfoTmp{3};
        a0(nf,nn) = cInfoTmp{4};
        del(nf,nn) = cInfoTmp{5};
        polya = cInfoTmp{6};
        polyp = cInfoTmp{7};
        calA0(nf,nn) = (nv(nf,nn)*l0(nf,nn))^2/(4*pi*a0(nf,nn));
        calA(nf,nn) = (polyp^2)/(4.0*pi*polya);

        % get vertex positions
        fline = fgetl(fid);     % goes to next line in file, skips header
        vPosTmp = textscan(fid,'VERTP %*f %f %f %f %f %f %f',NVERT); fline = fgetl(fid);     % goes to next line in file
        
        % parse data
        xposTmp = vPosTmp{1};
        yposTmp = vPosTmp{2};
        xvelTmp = vPosTmp{3};
        yvelTmp = vPosTmp{4};
        xfrcTmp = vPosTmp{5};
        yfrcTmp = vPosTmp{6};

        % save in cell
        xpos{nf,nn} = xposTmp;
        ypos{nf,nn} = yposTmp;
        xvel{nf,nn} = xvelTmp;
        yvel{nf,nn} = yvelTmp;
        xfrc{nf,nn} = xfrcTmp;
        yfrc{nf,nn} = yfrcTmp;
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
        NCELLStmp       = NCELLStmp{1};
        
        % read in packing fraction
        phi0            = textscan(fid,'PACKF %f',1);                   
        phi0            = phi0{1};
        fline = fgetl(fid);
        
        % update box size
        Ltmp            = textscan(fid,'BOXSZ %f %f',1);
        L(nf,1)         = Ltmp{1};
        L(nf,2)         = Ltmp{2};
        
         % skip next empty line
        fline = fgetl(fid);
        
        % get new contact info
        cijtmp      = textscan(fid,'NCONTS %f %f',1);
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    phi(nf:end) = [];
    Ncc(nf:end) = [];
    Nvv(nf:end) = [];
    xpos(nf:end,:) = [];
    ypos(nf:end,:) = [];
    xvel(nf:end,:) = [];
    yvel(nf:end,:) = [];
    xfrc(nf:end,:) = [];
    yfrc(nf:end,:) = [];
    nv(nf:end,:) = [];
    cxpos(nf:end,:) = [];
    cypos(nf:end,:) = [];
    l0(nf:end,:) = [];
    a0(nf:end,:) = [];
    del(nf:end,:) = [];
    calA(nf:end,:) = [];
    calA0(nf:end,:) = [];
    L(nf:end,:) = [];
end

% close position file
fclose(fid);

% store cell pos data into struct
cellJamData             = struct('NFRAMES',NFRAMES,'NCELLS',NCELLS);
cellJamData.phi         = phi;
cellJamData.Ncc         = Ncc;
cellJamData.Nvv         = Nvv;
cellJamData.xpos        = xpos;
cellJamData.ypos        = ypos;
cellJamData.xvel        = xvel;
cellJamData.yvel        = yvel;
cellJamData.xfrc        = xfrc;
cellJamData.yfrc        = yfrc;
cellJamData.nv          = nv;
cellJamData.cxpos       = cxpos;
cellJamData.cypos       = cypos;
cellJamData.l0          = l0;
cellJamData.a0          = a0;
cellJamData.del         = del;
cellJamData.calA        = calA;
cellJamData.calA0       = calA0;
cellJamData.L           = L;

end