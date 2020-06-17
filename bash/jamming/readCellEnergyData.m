function cellEnergyData = readCellEnergyData(fstr)
%% FUNCTION to read in cell energy data given file string

% print info to console
finfo = dir(fstr);
fprintf('Reading in %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% open file stream
fid = fopen(fstr);

% read in first line to discern number of particles
firstline = fgetl(fid);
firstdata = sscanf(firstline,'%f');

% parse first line of data

% non system size dependent
uintTmp = firstdata(1);
utotTmp = firstdata(2);
ktotTmp = firstdata(3);
frmsTmp = firstdata(4);
sxxTmp = firstdata(5);
sxyTmp = firstdata(6);
syyTmp = firstdata(7);
phiTmp = firstdata(8);
relaxItTmp = firstdata(9);

% system size dependent data
calAdata = firstdata(10:end);
sz = length(calAdata);
NCELLS = sz/2;
calA0inds = 1:2:sz-1;
calAinds = 2:2:sz;
calA0Tmp = calAdata(calA0inds);
calATmp = calAdata(calAinds);

% get rest of data
enfrmt = '%f %f %f %f %f %f %f %f %f';
calAfrmt = repmat(' %f %f',1,NCELLS);
endata = textscan(fid,[enfrmt calAfrmt]);
uintRest = endata{1};
utotRest = endata{2};
ktotRest = endata{3};
frmsRest = endata{4};
sxxRest = endata{5};
sxyRest = endata{6};
syyRest = endata{7};
phiRest = endata{8};
relaxItRest = endata{9};
NFRAMES = length(uintRest) + 1;

calA0Rest = endata(9 + calA0inds);
calARest = endata(9 + calAinds);

% initialize arrays
uint = zeros(NFRAMES,1);
utot = zeros(NFRAMES,1);
ktot = zeros(NFRAMES,1);
frms = zeros(NFRAMES,1);
stress = zeros(NFRAMES,3);
phi = zeros(NFRAMES,1);
relaxIt = zeros(NFRAMES,1);
calA = zeros(NFRAMES,NCELLS);
calA0 = zeros(NFRAMES,NCELLS);

% populate arrays
uint(1) = uintTmp;
uint(2:end) = uintRest;

utot(1) = utotTmp;
utot(2:end) = utotRest;

ktot(1) = ktotTmp;
ktot(2:end) = ktotRest;

frms(1) = frmsTmp;
frms(2:end) = frmsRest;

stress(1,1) = sxxTmp;
stress(2:end,1) = sxxRest;

stress(1,2) = sxyTmp;
stress(2:end,2) = sxyRest;

stress(1,3) = syyTmp;
stress(2:end,3) = syyRest;

phi(1) = phiTmp;
phi(2:end) = phiRest;

relaxIt(1) = relaxItTmp;
relaxIt(2:end) = relaxItRest;

for nn = 1:NCELLS
    calA(1,nn) = calATmp(nn);
    calA(2:end,nn) = calARest{nn};
    
    calA0(1,nn) = calA0Tmp(nn);
    calA0(2:end,nn) = calA0Rest{nn};
end


%% Make output into struct

cellEnergyData = struct('NFRAMES',NFRAMES);
cellEnergyData.NCELLS = NCELLS;
cellEnergyData.uint = uint;
cellEnergyData.utot = utot;
cellEnergyData.ktot = ktot;
cellEnergyData.frms = frms;
cellEnergyData.stress = stress;
cellEnergyData.phi = phi;
cellEnergyData.relaxIt = relaxIt;
cellEnergyData.calA = calA;
cellEnergyData.calA0 = calA0;


end