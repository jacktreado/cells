function gel2DSimPostprocess(varargin)
%% FUNCTION to process a gelation simulation after completion
% 
% 
% INPUTS:
% 
%   REQUIRED:
%   -- simulationFile: string that gives full path to *.pos file from simulation
%   -- saveFile: string that gives full path to .mat file with sim data
%       * Virial expression for stresses
%       * Force network over time
%       * RMS Forces/Energy over time
%       * packing fraction over time
% 
%   OPTIONAL:
%   -- simMovieFile: string that gives full path to .mp4 movie from
%       simulation (only run during gelation simulation)
% 
% NOTE ABOUT VIRIAL STRESSES: as of 09/27/2019, net force on each vertex
% is output in the .pos file, so in order to calculate the force between
% two cells we assume that each vertex is only interacting with vertices on
% 1 and only 1 other cell. THIS MAY NOT BE TRUE NEAR phi ~ 1 OR when a is
% large, so be aware that pairwise forces and therefore virial stresses are
% an approximation; to do them correctly, need to output pairwise forces
% between all cells during the entire simulation.
% 
% 
% ! ! CHANGE ABOVE COMMENT IF UPDATES TO CODE CHANGES THIS FACT ! !
% 

% check function arguments
if nargin < 2
    error('Need at least two arguments to function\n');
elseif nargin == 2
    % get arguments
    simulationFile  = varargin{1};
    saveFile        = varargin{2};
    writeToVideo    = 0;
    
    % test error in file names
    if ~strcmp(simulationFile(end-3:end),'.pos') && ~strcmp(simulationFile(end-4:end),'.test')
        error('Simulation file not *.pos or *.test, ending...\n');
    end
    if ~strcmp(saveFile(end-3:end),'.mat')
        error('Save file not *.mat, ending...\n');
    end
    
    % print if successful
    fprintf('Found 2 arguments to gelSimPostprocess function, so reading in from sim file %s and saving data to data file %s\n',simulationFile,saveFile);
elseif nargin == 3
    % get arguments
    simulationFile  = varargin{1};
    saveFile        = varargin{2};
    simMovieFile    = varargin{3};
    writeToVideo    = 1;
    
    % test error in file names
    if ~strcmp(simulationFile(end-3:end),'.pos') && ~strcmp(simulationFile(end-4:end),'.test')
        error('Simulation file not *.pos or *.test, ending...\n');
    end
    if ~strcmp(saveFile(end-3:end),'.mat')
        error('Save file not *.mat, ending...\n');
    end
    if ~strcmp(simMovieFile(end-3:end),'.mp4')
        error('Movie file not *.mp4, ending...\n');
    end
    
    % print if successful
    fprintf('Found 3 arguments to gelSimPostprocess function, so reading in from sim file %s and saving data to data file %s, \n\tand movie to movie file %s\n',simulationFile,saveFile,simMovieFile);
else
    error('4 or more arguments passed to gelSimPostprocess function, which is too many, ending.');
end
    
%% Load in data

% open file stream
fid = fopen(simulationFile);

% read in sim details
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);    NCELLS = NCELLS{1};
NFRAMES     = textscan(fid,'NUMFR %f',1);                   NFRAMES = NFRAMES{1};
BOXSZ       = textscan(fid,'BOXSZ %f',1);                     BOXSZ = BOXSZ{1};
fline       = fgetl(fid);

% initialize arrays to store per frame data
frameID = zeros(NFRAMES,1);
cxpos   = zeros(NFRAMES,NCELLS);
cypos   = zeros(NFRAMES,NCELLS);
l0      = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
calA    = zeros(NFRAMES,NCELLS);
xpos    = cell(NFRAMES,NCELLS);
ypos    = cell(NFRAMES,NCELLS);
xvel    = cell(NFRAMES,NCELLS);
yvel    = cell(NFRAMES,NCELLS);
xfrc    = cell(NFRAMES,NCELLS);
yfrc    = cell(NFRAMES,NCELLS);
nv      = zeros(NFRAMES,NCELLS);

% number of frames found
nf = 1;

while ~feof(fid)
    % start at NEWFR
    fline = fgetl(fid);
    newFrameStr = sscanf(fline,'%s %*s');
    if strcmp(newFrameStr,'STERM')
        fprintf('End of simulation found, NFRAMES = %d\n',nf);
        break;
    elseif ~strcmp(newFrameStr,'NEWFR') && ~strcmp(newFrameStr,'STERM')
        fprintf('File read-in not timed properly, does not start on new frame, newFrameStr = %s\n',newFrameStr);
        error('fix read structure');
    end
    
    % save frameID
    frameIDtmp = textscan(fid,'FRAME %f',1);
    frameID(nf) = frameIDtmp{1};

    % get number of vertices
    for nn = 1:NCELLS
        nvertTmp = textscan(fid,'NVERT %f',1); fline = fgetl(fid);     % goes to next line in file
        NVERT = nvertTmp{1};
        nv(nf,nn) = NVERT;
       
        % get cell pos and asphericity
        cInfoTmp = textscan(fid,'CELLP %f %f %f %f %f %f',1);   fline = fgetl(fid);     % goes to next line in file
        cxpos(nf,nn) = cInfoTmp{2};
        cypos(nf,nn) = cInfoTmp{3};
        l0(nf,nn) = cInfoTmp{4};
        a0(nf,nn) = cInfoTmp{5};
        calA(nf,nn) = cInfoTmp{6};

        % get vertex positions
        fline = fgetl(fid);     % goes to next line in file, skips header
        vPosTmp = textscan(fid,'VERTP %f %f %f %f %f %f %f',NVERT); fline = fgetl(fid);     % goes to next line in file
        
        % parse data
        xposTmp = vPosTmp{2};
        yposTmp = vPosTmp{3};
        xvelTmp = vPosTmp{4};
        yvelTmp = vPosTmp{5};
        xfrcTmp = vPosTmp{6};
        yfrcTmp = vPosTmp{7};

        % save in cell
        xpos{nf,nn} = xposTmp;
        ypos{nf,nn} = yposTmp;
        xvel{nf,nn} = xvelTmp;
        yvel{nf,nn} = yvelTmp;
        xfrc{nf,nn} = xfrcTmp;
        yfrc{nf,nn} = yfrcTmp;
    end
    
    % increment frame counter
    nf = nf + 1;
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    cxpos(nf:end,:) = [];
    cypos(nf:end,:) = [];
    xpos(nf:end,:) = [];
    ypos(nf:end,:) = [];
    xvel(nf:end,:) = [];
    yvel(nf:end,:) = [];
    xfrc(nf:end,:) = [];
    yfrc(nf:end,:) = [];
    l0(nf:end,:) = [];
    a0(nf:end,:) = [];
    calA(nf:end,:) = [];
    nv(nf:end,:) = [];
end

% close position file
fclose(fid);


%% Calculate pairwise forces on every frame

% gamma function to grab index of pairwise interaction
pw = @(i,j) (NCELLS*(i-1) + (j-1) - i*(i + 1)/2 + 1);

% number of pairwise forces
NPW = NCELLS*(NCELLS - 1)/2;

% initialize pairwise separation and force matrices
xpw = zeros(NFRAMES,NPW);
ypw = zeros(NFRAMES,NPW);
xfpw = zeros(NFRAMES,NPW);
yfpw = zeros(NFRAMES,NPW);

% initialize virial stresses in each frame
virialStress = zeros(NFRAMES,4);

% initialize root mean squared force in each frame
rmsForce = zeros(NFRAMES,1);

% cell to keep track of which vertices are interacting with what cell
NVERTS = nv(1,:);
vertexInt = cell(NCELLS,1);
for ci = 1:NCELLS
    vertexInt{ci} = zeros(NVERTS(ci),1);
end

% loop over frames, calculate pairwise forces
for ff = 1:NFRAMES
    % get vertex positions and forces
    cx  = cxpos(ff,:);
    cy  = cypos(ff,:);
    x   = xpos(ff,:);
    y   = ypos(ff,:);
    fx  = xfrc(ff,:);
    fy  = yfrc(ff,:);
    
    % loop over cells and verts, determine which cells they interact with
    for ci = 1:NCELLS
        % get vertex positions
        vx = x{ci};
        vy = y{ci};
        
        % vertex interaction vector
        vertIntTmp = vertexInt{ci};
        
        % loop over verts
        for vv = 1:NVERTS(ci)
            % vertex position
            vxTmp = vx(vv);
            vyTmp = vy(vv);
            
            % max dot product 
            cmax = 0.0;
            
            % get distance to other cells (ONLY NEED TO CHECK ci+1 AND UP
            % BECAUSE WHEN CALCULATING FORCE, VERTICES ON CELLS EARLIER IN
            % THE LIST WILL COMPENSATE FOR MISSING FORCE INFORMATION)
            for cj = ci+1:NCELLS
                % get center-to-center distance (J to I for dot product)
                cij = [cx(ci), cy(ci)] - [cx(cj), cy(cj)];
                cij = cij - BOXSZ*round(cij./BOXSZ);
                
                % store in separation vectors
                xpw(ff,pw(ci,cj)) = cij(1);
                ypw(ff,pw(ci,cj)) = cij(2);
                
                % get vector to other cell (J to I for dot product)
                rij =  [vxTmp, vyTmp] - [cx(cj), cy(cj)];
                rij = rij - BOXSZ*round(rij./BOXSZ);
                
                % get dot product between vertex- and cell-cell connectors
                if sum(cij.*rij) > cmax
                    vertIntTmp(vv) = cj;
                    cmax = sum(cij.*rij);
                end
            end
        end
        
        % store vertex interactions
        vertexInt{ci} = vertIntTmp;
    end
    
    % calculate pairwise forces
    for ci = 1:NCELLS
        % interaction information for cell ci
        vertIntTmp = vertexInt{ci};
        
        % forces
        vfx = fx{ci};
        vfy = fy{ci};
        
        % compute rms force
        for vv = 1:NVERTS(ci)
            rmsForce(ff) = rmsForce(ff) + vfx(vv)^2 + vfy(vv)^2;
        end
        
        % loop over pairs of cells
        for cj = ci+1:NCELLS
            % sum forces on ci due to cj
            cjInds = (vertIntTmp == cj);
            if sum(cjInds) == 0
                xfpw(ff,pw(ci,cj)) = 0.0;
                yfpw(ff,pw(ci,cj)) = 0.0;
            else
                xfpw(ff,pw(ci,cj)) = sum(vfx(cjInds));
                yfpw(ff,pw(ci,cj)) = sum(vfy(cjInds));
            end
            
            % add to virial stress
            virialStress(ff,1) = virialStress(ff,1) + 0.5*xpw(ff,pw(ci,cj))*xfpw(ff,pw(ci,cj));
            virialStress(ff,2) = virialStress(ff,2) + 0.5*ypw(ff,pw(ci,cj))*xfpw(ff,pw(ci,cj));
            virialStress(ff,3) = virialStress(ff,3) + 0.5*xpw(ff,pw(ci,cj))*yfpw(ff,pw(ci,cj));
            virialStress(ff,4) = virialStress(ff,4) + 0.5*ypw(ff,pw(ci,cj))*yfpw(ff,pw(ci,cj));
        end
    end
    
    % normalize rmsForce in this frame
    rmsForce(ff) = sqrt(rmsForce(ff))/(2*NCELLS);
end

% save time series data in matfile
save(saveFile,'rmsForce','xpw','ypw','xfpw','yfpw','cxpos','cypos','xpos','ypos','NVERTS','virialStress');

%% Record movie of simulation

if writeToVideo == 1
    % print message to console
    fprintf('Writing movie to file = %s\n',simMovieFile);
    
    % FRAME STEP
    NSTART = 1;
    NEND = NFRAMES;
    NSTEP = round(0.005*NFRAMES);
    if NSTEP < 1
        NSTEP = 1;
    end

    % figure window
    v = 1;

    % colors
    bidisperseScheme = jet(2);
    cellColor = zeros(NCELLS,3);
    cellColor(1:NCELLS/2,:) = repmat(bidisperseScheme(1,:),NCELLS/2,1);
    cellColor((NCELLS/2+1):end,:) = repmat(bidisperseScheme(2,:),NCELLS/2,1);
    lineColor = zeros(NCELLS,3);

    % write into movie
    vidObj = VideoWriter(simMovieFile,'MPEG-4');
    open(vidObj);

    figure(v), clf;

    fprintf('first visualizing on a loop...\n');
    for ff = NSTART:NSTEP:NEND
        % print to console
        fprintf('visualizing frame = %d, mean(a) = %f\n',frameID(ff),mean(calA(ff,:)));

        % visualize
        xprint = xpos(ff,:);
        yprint = ypos(ff,:);
        l0print = l0(ff,:);
        visualizeMultipleCells(v,NCELLS,BOXSZ,xprint,yprint,l0print,cellColor,lineColor);

        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end

    % close video object
    close(vidObj);
else
    fprintf('Skipping movie writing step, ending MATLAB function\n');
end


end