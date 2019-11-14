function qscompPostProcess(varargin)
%% FUNCTION to process a gelation simulation after completion
% 
% NOTE: this function will process data from multiple simulation seeds at a
% time
% 
% INPUTS:
% 
%   REQUIRED:
%   -- simulationLoc: folder with simulation data
%   -- simulationFile: string that gives simulation name
%       * Will then link to .pos, .en and .cm files
%   -- saveFile: string that gives full path to .mat file with sim data
%       * Virial expression for stresses
%       * Force network over time
%       * RMS Forces/Energy over time
%       * packing fraction over time
% 

% check function arguments
if nargin ~= 2
    fprintf('nargin = %d, but need nargin = 2, so ending.\n',nargin);
    error('Wrong number of arguments, ending\n');
else
    % get arguments
    simDir      = varargin{1};
    saveFile    = varargin{2};
    
    % check that simDir does not end in '/'
    if strcmp(simDir(end),'/')
        simDir = simDir(1:end-1);
    end
    
    % get name of sim directory in relative path
    [dataLoc,simDirName,~] = fileparts(simDir);
    
    % print if successful
    fprintf('Found 3 arguments to gelSimPostprocess function, so reading in from sim directory %s and saving data to data file %s\n',simDir,saveFile);
end

%% Get simulation directory info

% individual simulation folders
dirlist = dir([simDir '/qscomp*']);
NSIM = length(dirlist);
if NSIM == 0
    fprintf('No simulation subfolders found in sim dir %s\n',simDirName);
    error('No simulation subfolders found');
end

% == == == == == == == == ==
%   preliminary save data
% == == == == == == == == ==

% number of sims to neglect if need be
rmv         = false(NSIM,3);

% number of cells and frames
NCELLS      = zeros(NSIM,1);
NFRAMES     = zeros(NSIM,1);

% cell characteristics data
l0          = cell(NSIM,1);
a0          = cell(NSIM,1);
del         = cell(NSIM,1);
nv          = cell(NSIM,1);

% sim energy data
utot        = cell(NSIM,1);
uint        = cell(NSIM,1);
ktot        = cell(NSIM,1);
virial      = cell(NSIM,4);
phi         = cell(NSIM,1);

% sim contact network data
cm          = cell(NSIM,1);

% cell trajectory data
L           = cell(NSIM,1);
xpos        = cell(NSIM,1);
ypos        = cell(NSIM,1);
xvel        = cell(NSIM,1);
yvel        = cell(NSIM,1);
calA        = cell(NSIM,1);

%% Loop over simulation directories, parse data files and store data

for ss = 1:NSIM
    % print to console
    fprintf(' * * Processing data from %s\n',dirlist(ss).name);
    
    % specific simulation directory
    specificSimDir = [dirlist(ss).folder '/' dirlist(ss).name];
    
    % get file strings
    posfile = dir([specificSimDir '/*.pos']);
    enfile = dir([specificSimDir '/*.energy']);
    cmfile = dir([specificSimDir '/*.cm']);
    
    % check to see if position file exists
    if isempty(posfile) || posfile.bytes == 0
        fprintf('\t position file does not exist...\n')
        rmv(ss,:) = true(1,3);
        continue;
    else
        fprintf('\t position file %s does exists! parsing...\n',posfile.name);
        posfilestr = [posfile.folder '/' posfile.name];
    end
    
    % process position file
    [NFRAMES(ss), NCELLS(ss), L{ss}, xpos{ss}, ypos{ss}, xvel{ss}, yvel{ss}, virial{ss}, nv{ss}, l0{ss}, a0{ss}, del{ss}, calA{ss}] = parsePosFile(posfilestr);
    
    % check to see if energy file exists
    if isempty(enfile) || enfile.bytes == 0
        fprintf('\t energy file cannot be read in...\n')
        checkEn = 0;
        rmv(ss,2) = true;
    else
        fprintf('\t energy file %s does exists!\n',enfile.name);
        checkEn = 1;
        enfilestr = [enfile.folder '/' enfile.name];
    end
    
    % process energy file
    if checkEn == 1
        fprintf('\t Parsing energy file\n');
        [utot{ss}, uint{ss}, ktot{ss}, phi{ss}] = parseEnFile(enfilestr);
    end
    
    if isempty(cmfile) || cmfile.bytes == 0
        fprintf('\t contact matrix file does not exist...\n')
        checkCM = 0;
        rmv(ss,3) = true;
    else
        fprintf('\t contact matrix file %s does exists!\n',cmfile.name);
        checkCM = 1;
        cmfilestr = [cmfile.folder '/' cmfile.name];
    end
    
     % process contact matrix file
    if checkCM == 1
        fprintf('\t Parsing contact matrix file\n');
        cm{ss} = parseCMFile(cmfilestr,NCELLS(ss));
    end
end

% save to save file
save([dataLoc '/' saveFile],'NFRAMES','NCELLS','L','xpos','ypos','xvel','yvel',...
    'virial','nv','l0','a0','del','calA',...
    'utot','uint','ktot','phi',...
    'cm','rmv','NSIM');


end






%% FUNCTION TO PARSE INDIVIDUAL .pos FILE AND OUTPUT SIM VALUES


function [NFRAMES, NCELLS, L, xpos, ypos, xvel, yvel, virial, nv, l0, a0, del, calA] = parsePosFile(filestr)
% get file info
finfo = dir(filestr);
fname = finfo.name;
fsize = finfo.bytes;

% check extension
[~, ~, ext] = fileparts(filestr);

% print file info
if ~strcmp(ext,'.pos')
    fprintf('ERROR: in parsePosFile function, FILE %s DOES NOT HAVE CORRECT EXTENSION, ENDING\n',fname);
    error('input file has incorrect extension');
else
    fprintf('Processing filestr = %s, is of size %.3f MB\n',fname,fsize/1e6);
end

% open file stream
fid = fopen(filestr);

% read in sim details from first frame
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);   
NCELLS      = NCELLS{1};
Ltmp        = textscan(fid,'BOXSZ %f %f',1);    
fline       = fgetl(fid);

% variables to store data
NFRAMES = 5e4;
xpos    = cell(NFRAMES,NCELLS);
ypos    = cell(NFRAMES,NCELLS);
xvel    = cell(NFRAMES,NCELLS);
yvel    = cell(NFRAMES,NCELLS);
virial  = zeros(NFRAMES,4);
nv      = zeros(NCELLS,1);
l0      = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
del     = zeros(NFRAMES,NCELLS);
calA    = zeros(NFRAMES,NCELLS);
L       = zeros(NFRAMES,2);

% number of frames found
nf = 1;

% loop over lines of file, parse frame by frame
while ~feof(fid)
    % get virial stress information
    vstress = textscan(fid,'VRIAL %f %f %f %f',1);
    virial(nf,1) = vstress{1};
    virial(nf,2) = vstress{2};
    virial(nf,3) = vstress{3};
    virial(nf,4) = vstress{4};
    
    % get box length
    L(nf,1) = Ltmp{1};
    L(nf,2) = Ltmp{2};
    
    % discard header line
    fline = fgetl(fid);
    
    % get info about deformable particle
    for nn = 1:NCELLS
        nvertTmp = textscan(fid,'NVERT %f',1); fline = fgetl(fid);     % goes to next line in file
        NVERT = nvertTmp{1};
        if nf == 1
            nv(nn) = NVERT;
        end
       
        % get cell pos and asphericity
        cInfoTmp    = textscan(fid,'CELLP %*f %f %f %f %f %f %f',1);   fline = fgetl(fid);     % goes to next line in file
        cxpos       = cInfoTmp{1};
        cypos       = cInfoTmp{2};
        l0(nf,nn)   = cInfoTmp{3};
        a0(nf,nn)   = cInfoTmp{4};
        del(nf,nn)  = cInfoTmp{5};
        calA(nf,nn) = cInfoTmp{6};

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
    end
    
    % increment frame count
    nf = nf + 1;
    
    % read in/discard trash information
    % NOTE: IF BOX GEOMETRY/NUMBER OF PARTICLES/ANYTHING ELSE CHANGES
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
        NCELLStmp       = textscan(fid,'NUMCL %f',1,'HeaderLines',1);
        
        % NOTE: IF BOX SIZE CHANGES, UPDATE L HERE
        Ltmp            = textscan(fid,'BOXSZ %f %f',1);
        L(nf,1)         = Ltmp{1};
        L(nf,2)         = Ltmp{2};
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    xpos(nf:end,:) = [];
    ypos(nf:end,:) = [];
    xvel(nf:end,:) = [];
    yvel(nf:end,:) = [];
    virial(nf:end,:) = [];
    nv(nf:end,:) = [];
    l0(nf:end,:) = [];
    a0(nf:end,:) = [];
    del(nf:end,:) = [];
    calA(nf:end,:) = [];
    L(nf:end,:) = [];
end

fprintf('finished parsing all %d frames\n',NFRAMES);

% close file stream
fclose(fid);

end







%% FUNCTION TO PARSE INDIVIDUAL .en FILE AND OUTPUT ENERGY VALUES

function [utot, uint, ktot, phi] = parseEnFile(enfilestr)

% get file info
finfo = dir(enfilestr);
fname = finfo.name;
fsize = finfo.bytes;

% check extension
[~, ~, ext] = fileparts(enfilestr);

% print file info
if ~strcmp(ext,'.energy')
    fprintf('ERROR: in parseEnFile function, FILE %s DOES NOT HAVE CORRECT EXTENSION, ENDING\n',fname);
    error('input file has incorrect extension');
else
    fprintf('Processing filestr = %s, is of size %.3f MB\n',fname,fsize/1e6);
end

% open file stream
fid = fopen(enfilestr);

% read
edata   = textscan(fid,'%f %f %f %*f %*f %*f %*f %f');
uint    = edata{1};
utot    = edata{2};
ktot    = edata{3};
phi     = edata{4};

% close file stream
fclose(fid);

end





%% FUNCTION TO PARSE INDIVIDUAL .cm FILE AND OUTPUT ENERGY VALUES

function cm = parseCMFile(cmfilestr,NCELLS)

% get file info
finfo = dir(cmfilestr);
fname = finfo.name;
fsize = finfo.bytes;

% check extension
[~, ~, ext] = fileparts(cmfilestr);

% print file info
if ~strcmp(ext,'.cm')
    fprintf('ERROR: in parseCMFile function, FILE %s DOES NOT HAVE CORRECT EXTENSION, ENDING\n',fname);
    error('input file has incorrect extension');
else
    fprintf('Processing filestr = %s, is of size %.3f MB\n',fname,fsize/1e6);
end

% open file stream
fid = fopen(cmfilestr);

% contact matrix
NCONTACTS = NCELLS*(NCELLS-1)/2;
NFRAMES = 1e6;
cm = zeros(NFRAMES,NCONTACTS);

% number of frames read
nf = 1;

% loop over file
while ~feof(fid)
    % get line as string
    fline = fgetl(fid);

    % parse string for contact vals
    cm(nf,:) = sscanf(fline,'%f');

    % increment frame counter
    nf = nf + 1;
end

% delete extra entries
cm(nf:end,:) = [];

% close file stream
fclose(fid);

end








