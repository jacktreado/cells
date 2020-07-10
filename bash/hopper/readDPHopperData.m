function dpHopperPosData = readDPHopperData(fstr)
%% function to read in position data from clogging simulation in hopper
% geometry

% print info to console
finfo = dir(fstr);
fprintf('Reading in clogging pos data: %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% open file stream
fid = fopen(fstr);

% read in sim details from first frame
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);   
NCELLS      = NCELLS{1};

% get number of interaction pairs
NPAIRS = 0.5*NCELLS*(NCELLS-1);
cmfrmt = repmat('%f ',1,NPAIRS);

hopperData  = textscan(fid,'HOPPR %f %f %f',1);
w0          = hopperData{1};
w           = hopperData{2};
th          = hopperData{3};

% initialize arrays to store per frame data
NFRAMES = 5e4;
wallfrc = zeros(NFRAMES,4);
Ncc     = zeros(NFRAMES,1);
Nvv     = zeros(NFRAMES,1);
cm      = zeros(NFRAMES,NPAIRS);
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

% number of frames found
nf = 1;

% loop over lines of file, parse frame by frame
while ~feof(fid)
    % get virial stress information
    vstress         = textscan(fid,'WLFRC %f %f %f %f',1);
    wallfrc(nf,1)   = vstress{1};
    wallfrc(nf,2)   = vstress{2};
    wallfrc(nf,3)   = vstress{3};
    wallfrc(nf,4)   = vstress{4};
    
    % get contact information
    cntinfo         = textscan(fid,'NCNTS %f %f',1);
    Ncc(nf)         = cntinfo{1};
    Nvv(nf)         = cntinfo{2};
    
    % get contact matrix
    cmdata          = textscan(fid,['CMMAT' cmfrmt],1);
    cm(nf,:)        = cell2mat(cmdata);
    
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
        
        % read in hopper data for next frame
        hopperData  = textscan(fid,'HOPPR %f %f %f',1);
        w0tmp       = hopperData{1};
        wtmp        = hopperData{2};
        thtmp       = hopperData{3};
        
         % skip next empty line
        fline = fgetl(fid);
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    wallfrc(nf:end,:) = [];
    Ncc(nf:end) = [];
    Nvv(nf:end) = [];
    cm(nf:end,:) = [];
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
end

% close position file
fclose(fid);

% calculate length of hopper nozzle
sb = l0(1,1)*del(1,1);
L = 0.5*(w0 - w)*tan(th) + 0.5*sb*((1/cos(th)) + 1 - tan(th));

% store cell pos data into struct
dpHopperPosData             = struct('NFRAMES',NFRAMES,'NCELLS',NCELLS);
dpHopperPosData.w           = w;
dpHopperPosData.w0          = w0;
dpHopperPosData.th          = th;
dpHopperPosData.L           = L;
dpHopperPosData.wallfrc     = wallfrc;
dpHopperPosData.Ncc         = Ncc;
dpHopperPosData.Nvv         = Nvv;
dpHopperPosData.cm          = cm;
dpHopperPosData.xpos        = xpos;
dpHopperPosData.ypos        = ypos;
dpHopperPosData.xvel        = xvel;
dpHopperPosData.yvel        = yvel;
dpHopperPosData.xfrc        = xfrc;
dpHopperPosData.yfrc        = yfrc;
dpHopperPosData.nv          = nv;
dpHopperPosData.cxpos       = cxpos;
dpHopperPosData.cypos       = cypos;
dpHopperPosData.l0          = l0;
dpHopperPosData.a0          = a0;
dpHopperPosData.del         = del;
dpHopperPosData.calA        = calA;
dpHopperPosData.calA0       = calA0;
dpHopperPosData.L           = L;


end