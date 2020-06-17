function cellVDOSData = readCellVDOS(vdosfstr,NFRAMES)
%% Function to read in VDOS data from cell file
% NOTE: assume NFRAMES known externally, in future
% add varargin to use eof to read in when NFRAMES not known

% open file stream
fid = fopen(vdosfstr);

finfo = dir(vdosfstr);
fprintf('Reading in %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% initialize output data
NDOF            = zeros(NFRAMES,1);
shapeStEvals    = cell(NFRAMES,1);
contactStEvals  = cell(NFRAMES,1);
totalStEvals    = cell(NFRAMES,1);
totalEvals      = cell(NFRAMES,1);
totalEvecs      = cell(NFRAMES,1);

% loop over frames, grab vdos info
for ff = 1:NFRAMES
    % scan in data
    data = textscan(fid,'%d',1);
    NDOF(ff) = data{1};

    % read in eigenvalues of stiffness matrix
    data = textscan(fid,'%f',NDOF(ff));
    shapeStEvals{ff} = data{1};
    
    % read in eigenvalues of stress matrix
    data = textscan(fid,'%f',NDOF(ff));
    contactStEvals{ff} = data{1};
    
    % read in eigenvalues of total matrix
    data = textscan(fid,'%f',NDOF(ff));
    totalStEvals{ff} = data{1};
    
    % read in eigenvalues of total matrix
    data = textscan(fid,'%f',NDOF(ff));
    totalEvals{ff} = data{1};

    % read in eigenvectors
    frmt = repmat('%f ',1,NDOF(ff));
    data = textscan(fid,frmt,NDOF(ff));
    evecstmp = zeros(NDOF(ff));
    for vv = 1:NDOF(ff)
        evecstmp(:,vv) = data{vv};
    end
    totalEvecs{ff} = evecstmp;
    
%     % read in energy due to perturbations
%     data = textscan(fid,'%d',1);
%     NSTEPS(ff) = data{1}; 
%     
%     data = textscan(fid,'%f',1);
%     U0(ff) = data{1};
%     
%     frmt = repmat('%f ',1,NSTEPS(ff));
%     data = textscan(fid,frmt,1);
%     detmp = zeros(1,NSTEPS(ff));
%     for ss = 1:NSTEPS(ff)
%         detmp(ss) = data{ss};
%     end
%     de{ff} = detmp;
%     
%     data = textscan(fid,frmt,NDOF(ff));
%     dUtmp = zeros(NDOF(ff),NSTEPS(ff));
%     for vv = 1:NSTEPS(ff)
%         dUtmp(:,vv) = data{vv};
%     end
%     Ut{ff} = dUtmp;
end

% close file object
fclose(fid);

% save data
cellVDOSData = struct('NDOF',NDOF);
cellVDOSData.shapeStEvals = shapeStEvals;
cellVDOSData.contactStEvals = contactStEvals;
cellVDOSData.totalStEvals = totalStEvals;
cellVDOSData.totalEvals = totalEvals;
cellVDOSData.totalEvecs = totalEvecs;

end