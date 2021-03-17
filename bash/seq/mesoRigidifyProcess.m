function mesoRigidifyProcess(fpattern,savefstr)
%% FUNCTION to process meso rigidify simulations

% get list of files
fprintf('Reading files from %s\n',fpattern);
flist = dir([fpattern '*.pos']);
NF = length(flist);
if NF == 0
    error('No files found.');
else
    fprintf('Found files NF = %d\n',NF);
end

%% Loop over files

% -- bin face info

% phi0 bins
pbins       = 20;
phi0BE      = linspace(0.95,0.385,pbins+1);
phi0BC      = 0.5*(phi0BE(2:end) + phi0BE(1:end-1));

% face bins
fbins       = 30;
faceBC      = (1:fbins)' + 2;

% face counts
faceCounts  = zeros(fbins,pbins);

% -- save simulation info
LList       = zeros(NF,2);
nvList      = cell(NF,1);
phi0List    = cell(NF,1);
calAList    = cell(NF,1);
calA0List   = cell(NF,1);
cijList     = cell(NF,1);
NFvalsList  = cell(NF,1);

for ii = 1:NF
    % get pos file str
    fstr = [flist(ii).folder '/' flist(ii).name];
    
    % contact network information
    ctcstr = [fstr(1:end-4) '.ctc'];
    
    % read in data
    mesophyllTrajectoryData = readMesophyllData(fstr);

    % get number of frames
    NFRAMES     = mesophyllTrajectoryData.NFRAMES;
    NCELLS      = mesophyllTrajectoryData.NCELLS;
    nv          = mesophyllTrajectoryData.nv(1,:);
    L           = mesophyllTrajectoryData.L(1,:);
    Lx          = L(1);
    Ly          = L(2);

    % load in ctc data
    cijInfo     = cell(NFRAMES,1);
    NCTCS       = 0.5*NCELLS*(NCELLS-1);
    ltInds      = find(tril(ones(NCELLS),-1));
    frmt        = repmat('%f ',1,NCTCS);
    fid         = fopen(ctcstr);
    for ff = 1:NFRAMES
        cijtmp          = zeros(NCELLS);
        ctmp            = textscan(fid,frmt,1);
        ctmp            = cell2mat(ctmp);
        cijtmp(ltInds)  = ctmp;
        cijtmp          = cijtmp + cijtmp';
        cijInfo{ff}     = cijtmp;
    end
    fclose(fid);
    
    % Loop over frames, get phi0, calA0
    l0mat = mesophyllTrajectoryData.l0;
    a0mat = mesophyllTrajectoryData.a0;
    calA0 = mesophyllTrajectoryData.calA0;
    phi0 = zeros(NFRAMES,1);
    for ff = 1:NFRAMES
        phi0(ff) = sum(a0mat(ff,:) + (0.25*(l0mat(ff,:).^2).*(0.5*nv - 1)))/(Lx*Ly);
    end
    
    % save data
    LList(ii,:)     = L;
    phi0List{ii}    = phi0;
    calA0List{ii}   = calA0;
    nvList{ii}      = nv;
    cijList{ii}     = cijInfo;
    
    % periodic images info
    NP              = 9*NCELLS;
    
    % Loop over frames, get particle shapes and void polygons
    calA = zeros(NFRAMES,NCELLS);
    NFvals = zeros(NFRAMES,2);
    for ff = 1:NFRAMES
        % get contacts
        cij = cijInfo{ff};
        
        % get index of phi0 for below
        phi0BinIdx = (phi0(ff) > phi0BE(1:end-1) & phi0(ff) < phi0BE(2:end));
    
        % cell positions at frame ff
        xpos = mesophyllTrajectoryData.xpos(ff,:);
        ypos = mesophyllTrajectoryData.ypos(ff,:);
        
        % loop over particles, compute instantaneous calA
        for nn = 1:NCELLS
            % positions
            xtmp = xpos{nn};
            ytmp = ypos{nn};
            
            % area
            atmp = polyarea(xtmp,ytmp);
            
            % indexing for perimeter
            ip1 = [2:nv(nn) 1];
            
            % perimeter
            lx = xtmp(ip1) - xtmp;
            ly = ytmp(ip1) - ytmp;
            l = sqrt(lx.^2 + ly.^2);
            ptmp = sum(l);
            
            % save instantaneous shape parameter
            calA(ff,nn) = ptmp^2/(4.0*pi*atmp);
        end
        
        % print to console
        fprintf('\t - - Computing for frame ff = %d/%d',ff,NFRAMES);
        
        % loop over particle centers
        fprintf('... image pos');
        P = zeros(NP,2);
        for nn = 1:NCELLS
            cx = mean(xpos{nn});
            cy = mean(ypos{nn});
            kbox = 1;
            for xx = -1:1
                for yy = -1:1
                    P((kbox-1)*NCELLS + nn,1) = cx + xx*Lx;
                    P((kbox-1)*NCELLS + nn,2) = cy + yy*Ly;
                    kbox = kbox + 1;
                end
            end
        end
        
        fprintf('... contacts');
        A = zeros(NP);
        for nn = 1:NCELLS
            cx = mean(xpos{nn});
            cy = mean(ypos{nn});
            for mm = nn+1:NCELLS
                if cij(nn,mm) > 0
                    xm = mean(xpos{mm});
                    ym = mean(ypos{mm});

                    dx = xm - cx;
                    imx = round(dx/Lx);
                    dx = dx - Lx*imx;

                    dy = ym - cy;
                    imy = round(dy/Ly);
                    dy = dy - Ly*imy;

                    for xx = -1:1
                        xi = 2 + xx;
                        for yy = -1:1
                            yi = 2 + yy;
                            bi = 3*(xi-1) + yi;

                            % add to periodic planar graph if not on exterior boundary
                            ytop = (imy == -1 && yy == 1);
                            xright = (imx == -1 && xx == 1);
                            ybottom = (imy == 1 && yy == -1);
                            xleft = (imx == 1 && xx == -1);

                            if (~ytop && ~ybottom && ~xleft && ~xright)
                                bimm = 3*(xi - imx - 1) + yi - imy;

                                nnInd = (bi-1)*NCELLS + nn;
                                mmInd = (bimm-1)*NCELLS + mm;
                                A(nnInd,mmInd) = 1;
                                A(mmInd,nnInd) = 1;
                            end
                        end
                    end
                end
            end
        end
        
        % get spatial graphs
        fprintf('... spatial graph\n');
        G = graph(A);
        obj = spatialgraph2D(G,P(:,1),P(:,2));
        pgon = polyshape(obj);
        
        vpos = {pgon.Vertices}';
        NFACES = size(vpos,1);
        NFMAIN = 0;
        for jj = 1:NFACES
            vptmp = vpos{jj};
            FV = length(vptmp);
            cx = mean(vptmp(:,1));
            cy = mean(vptmp(:,2));
            if cx > 0 && cx < Lx && cy > 0 && cy < Ly
                % increment number of faces in main
                NFMAIN = NFMAIN + 1;
                
                % add face info to face count hist
                faceCounts(FV-2,phi0BinIdx) = faceCounts(FV-2,phi0BinIdx) + 1;
            end
        end
        
        % print info about NFaces vs Euler characteristic
        nc = 0.5*sum(cij>0,'all');
        NFEULER = nc - NCELLS;
        fprintf('\t \t ** NFMAIN = %d, NFEULER = %d\n',NFMAIN,NFEULER);
        NFvals(ff,1) = NFEULER;
        NFvals(ff,2) = NFMAIN;
    end
    
    % save [article info
    calAList{ii} = calA;
    NFvalsList{ii} = NFvals;
end

%% Save to save string

save(savefstr,'flist','pbins','fbins','faceBC','phi0BC',...
    'LList','nvList','phi0List','calAList','calA0List',...
    'cijList','NFvalsList');

end