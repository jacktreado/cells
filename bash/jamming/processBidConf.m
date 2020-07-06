function processBidConf(simStr,saveStr)
%% FUNCTION to process compress-to-confluence simulations, extract VDOS, contacts, etc

% get list of files
flist = dir([simStr '*.vdos']);
NSIM = length(flist);
if NSIM == 0
    fprintf('No files found using str %s\n',simStr);
    error('No files found, ending.');
end

%% loop over simulations, extract data, save

% arrays
simList         = cell(NSIM,1);         % simulation file name list
NCELLSList      = zeros(NSIM,1);        % number of cells in each sim
NDOFList        = zeros(NSIM,1);        % total # of degrees of freedom in each sim
NFRAMEList      = zeros(NSIM,1);        % total # of frames for each sim
NvList          = cell(NSIM,1);         % # of vertices on each particle
LList           = zeros(NSIM,2);        % box lengths
a0List          = cell(NSIM,1);         % particle preferred area list
l0List          = cell(NSIM,1);         % vertex size list

% containers for extracted VDOS data
pList           = cell(NSIM,1);         % list of positive pressures
NvvList         = cell(NSIM,1);         % # of vertex-vertex contacts
NccList         = cell(NSIM,1);         % # of cell-cell contacts
NqList          = cell(NSIM,1);         % # of quartic modes, counted by stiffness matrix
lambdaMatList   = cell(NSIM,1);         % matrix of lambda values as a function of pressure
stiffEVMatList  = cell(NSIM,1);         % matrix of eigenvalues from stiffness matrix
vertexPrList    = cell(NSIM,1);         % vertex participation ratio as a function of pressure
cellPrList      = cell(NSIM,1);         % cell participation ratio as a function of pressure     
projList        = cell(NSIM,1);         % mode projection as a function of pressure (Dong's 3 directions, Frenet-Serret)

% containers for extracted energetic data
jFrameList      = cell(NSIM,1);         % map from energetic frame list to vdos frame list
enPList         = cell(NSIM,1);         % pressure from energetic data
enPhiList       = cell(NSIM,1);         % phi from energetic data
enCalAList      = cell(NSIM,2);         % calA from energetic data

% simulation slots to delete
simSkip         = false(NSIM,1);        

% loop over simulations
for ss = 1:NSIM
    % file information
    dirname                 = flist(ss).folder;
    vdosfname               = flist(ss).name;
    jamfname                = [vdosfname(1:end-5) '.pos'];
    enfname                 = [vdosfname(1:end-5) '.en'];
    
    jamfstr                 = [dirname '/' jamfname];
    enfstr                  = [dirname '/' enfname];
    vdosfstr                = [dirname '/' vdosfname];
    
    % test to make sure that all files exist and have data
    if exist(jamfstr,'file') > 0
        jamfinfo = dir(jamfstr);
        if jamfinfo.bytes == 0
            fprintf('jam file %s exists but is empty, skipping sim...\n',jamfstr);
            simSkip(ss) = true;
            continue;
        else
            if exist(enfstr,'file') > 0
                enfinfo = dir(enfstr);
                if enfinfo.bytes == 0
                    fprintf('energy file %s exists but is empty, skipping sim...\n',enfstr);
                    simSkip(ss) = true;
                    continue;
                else
                    if exist(vdosfstr,'file') > 0
                        vdosfinfo = dir(vdosfstr);
                        if vdosfinfo.bytes == 0
                            fprintf('vdos file %s exists but is empty, skipping sim...\n',vdosfstr);
                            simSkip(ss) = true;
                            continue;
                        else
                            simList{ss} = vdosfstr;
                        end
                    else
                        fprintf('vdos file %s does not exist in location, skipping sim...\n',vdosfstr);
                        simSkip(ss) = true;
                        continue;
                    end
                end
            else
                fprintf('energy file %s does not exist in location, skipping sim...\n',enfstr);
                simSkip(ss) = true;
                continue;
            end
        end
    else
        fprintf('jam file %s does not exist in location, skipping sim...\n',jamfstr);
        simSkip(ss) = true;
        continue;
    end
    
    % print to console
    fprintf('\t ** Reading in simulation data from %s...\n',vdosfname);
    
    % load data
    cellJamData             = readCellJamData(jamfstr);

    % parse time constant information
    NFRAMES                 = cellJamData.NFRAMES;
    NCELLS                  = cellJamData.NCELLS;
    L                       = cellJamData.L;
    nv                      = cellJamData.nv(1,:);
        
    % total number of degrees of freedom
    NVTOT                   = sum(nv);
    NDOF                    = 2*NVTOT;

    % set box size
    Lx                      = L(1,1);
    Ly                      = L(1,2);
    
    % store system data for this sim
    NCELLSList(ss)          = NCELLS;
    NDOFList(ss)            = NDOF;
    NvList{ss}              = nv;
    LList(ss,:)             = [Lx, Ly];

    % vdos data
    cellVDOSData            = readCellVDOS(vdosfstr,NFRAMES);

    % energy data
    cellEnergyData          = readCellEnergyData(enfstr);
    NENFRAMES               = cellEnergyData.NFRAMES;

    % determine which frames to skip
    if NENFRAMES ~= NFRAMES
        % align frames by phi
        phiEn = cellEnergyData.phi;
        phiC = cellJamData.phi;

        % assuming more energy frames than jamming frames
        jframes = false(NENFRAMES,1);
        rmvphi = false(NFRAMES,1);
        for ff = 1:NFRAMES
            [dmin,matchInd] = min(abs(phiC(ff) - phiEn));
            if ~jframes(matchInd)
                jframes(matchInd) = true;
            else
                if dmin > 1e-10
                    error('tried to double align a large difference, ending.');
                else
                    rmvphi(ff) = true;
                end
            end
        end
        phiC(rmvphi) = [];
        NFRAMES = length(phiC);

        if (sum(jframes) ~= NFRAMES)
            error('Frame alignment incorrect, ending.');
        end
    else
        jframes = true(NFRAMES,1);
    end
    
    fprintf('\t ** Parsing simulation data...\n');
    
    % get pressure data over frames, sort frames by pressure
    PCUT                    = 1e-10;
    stress                  = cellEnergyData.stress;
    Pinput                  = 0.5*(stress(jframes,1) + stress(jframes,3))/(NCELLS*Lx*Ly);
    [P,inds]                = sort(Pinput);
    ssPos                   = inds(P > PCUT);
    
    % parse and sort vdos data
    NvvInput                = cellJamData.Nvv;
    NccInput                = cellJamData.Ncc;
    xposInput               = cellJamData.xpos;
    yposInput               = cellJamData.ypos;
    l0Input                 = cellJamData.l0;
    a0Input                 = cellJamData.a0;
    lambdaInput             = cellVDOSData.totalEvals;
    evecsInput              = cellVDOSData.totalEvecs;
    totalStEvalsInput       = cellVDOSData.totalStEvals;

    NvvPos                  = NvvInput(ssPos);
    NccPos                  = NccInput(ssPos);
    lambdaPos               = lambdaInput(ssPos);
    stiffEPos               = totalStEvalsInput(ssPos);
    evecsPos                = evecsInput(ssPos);
    Ppos                    = P(P > PCUT);
    NPOSFRAMES              = sum(P > PCUT);
    xpos                    = xposInput(ssPos,:);
    ypos                    = yposInput(ssPos,:);
    l0Pos                   = l0Input(ssPos,:);
    a0Pos                   = a0Input(ssPos,:);
    
    NFRAMEList(ss)          = NPOSFRAMES;
    
    fprintf('\t ** Storing simulation data...\n');
    
    % store extracted data
    pList{ss}               = Ppos;
    NvvList{ss}             = NvvPos;
    NccList{ss}             = NccPos;
    l0List{ss}              = l0Pos;
    a0List{ss}              = a0Pos;
    
    fprintf('\t ** -- Eigenvalue matrix...\n');
    
    % make matrix of eigenvalues
    lambdaMatrix = zeros(NDOF,NPOSFRAMES);
    for pp = 1:NPOSFRAMES
        lambdaMatrix(:,pp) = lambdaPos{pp};
    end
    
    % make matrix of eigenvalues from stiffness matrix
    stiffEVMatrix = zeros(NDOF,NPOSFRAMES);
    for pp = 1:NPOSFRAMES
        stiffEVMatrix(:,pp) = stiffEPos{pp};
    end
    
    fprintf('\t ** -- Quartic modes...\n');
    
    % compute number of quartic modes from stiffness matrix
    nq = zeros(NPOSFRAMES,1);
    for pp = 1:NPOSFRAMES
        nq(pp) = sum(stiffEPos{pp} < 1e-8) - 2;
    end
    
    fprintf('\t ** -- Participation ratios...\n');
    
    % compute participation ratio at each pressure
    vertexPR = zeros(NDOF,NPOSFRAMES);
    cellPR = zeros(NDOF,NPOSFRAMES);
    for pp = 1:NPOSFRAMES

        % get eigenvector
        evecs = evecsPos{pp};

        % loop over modes (vertex participation ratio)
        for ww = 1:NDOF
            % get eigenvector
            evtmp = evecs(:,ww);

            % compute numerator for vertex pr
            prnum = sum(evtmp.^2)^2;

            % compute denominator for vertex pr
            prdenom = 0.0;
            for ii = 1:2:(NDOF-1)
                prdenom = prdenom + (evtmp(ii)^2 + evtmp(ii+1)^2)^2;
            end
            
            % compute numerator for cell pr
            cellprnum = 0.0;
            cellprdenom = 0.0;
            i0 = 0;
            for nn = 1:NCELLS
                xi = (i0 + 1):2:(i0 + nv(nn) - 1);
                yi = (i0 + 2):2:(i0 + nv(nn));
                i0 = nv(nn);
                
                cdx = sum(evtmp(xi))/nv(nn);
                cdy = sum(evtmp(yi))/nv(nn);
                
                cellprnum = cellprnum + (cdx^2 + cdy^2);
                cellprdenom = cellprdenom + (cdx^2 + cdy^2)^2; 
            end

            % compute vertex participation ratio
            vertexPR(ww,pp) = prnum/(NVTOT*prdenom);
            
            % compute cell participation ratio
            cellPR(ww,pp) = (cellprnum^2)/(NCELLS*cellprdenom);
        end
    end
    
    fprintf('\t ** -- Mode projections...\n');
    
    % compute projections for each pressure
    projections = cell(NPOSFRAMES,5);
    for pp = 1:NPOSFRAMES
        % swap eigenvector order for use in Dong's code
        evectmp = evecsPos{pp};
        evec = evectmp;
        
        % my indexing
        x0inds = 1:2:(NDOF-1);
        y0inds = 2:2:NDOF;

        % Dong's indexing
        xinds = 1:NVTOT;
        yinds = (NVTOT + 1):NDOF;

        evec(xinds,:) = evectmp(x0inds,:);
        evec(yinds,:) = evectmp(y0inds,:);

        % make total list of vertex degrees of freedom
        xall = zeros(NVTOT,1);
        yall = zeros(NVTOT,1);
        Dc = zeros(NCELLS,1);

        xtmp = xpos(pp,:);
        ytmp = ypos(pp,:);
        l0tmp = l0Pos(pp,:);
        last = 1;

        for nn = 1:NCELLS
            next = last + nv(nn) - 1;
            xall(last:next) = xtmp{nn};
            yall(last:next) = ytmp{nn};
            last = next + 1;

            % compute effective cell diameter
            Dc(nn) = l0tmp(nn)/sin(pi/nv(nn));
        end
        
        % use Dong's code
        V2_norm = ModeProj_DPM(NCELLS, nv, Dc, xall, yall, evec);
        
        % save projections for this pressure
        projections{pp,1} = V2_norm(1,:)';
        projections{pp,2} = V2_norm(2,:)';
        projections{pp,3} = V2_norm(3,:)';
        
        % Project eigenvectors onto frenet serret
        
        % compute FS directions
        tv = zeros(NDOF,1);
        sv = zeros(NDOF,1);
        last = 1;
        for nn = 1:NCELLS
            % vertex indices
            ip1 = [(2:nv(nn))'; 1];
            im1 = [nv(nn); (1:nv(nn)-1)'];
            
            % vertex coordinates for cell nn
            vx = xtmp{nn};
            vy = ytmp{nn};

            % directions for each vertex, normalized for each cell
            vl = [vx(ip1) - vx, vy(ip1) - vy];
            l = sqrt(sum(vl.^2,2));
            ul = vl./l;
            ti = ul(im1,:) + ul;
            ti = ti./sqrt(sum(ti.^2,2));
            si = ul(im1,:) - ul;
            si = si./sqrt(sum(si.^2,2));
            
            % store in array
            next = last + 2*nv(nn) - 1;
            xindnn = last:2:(next-1);
            yindnn = (last+1):2:next;
            tv(xindnn) = ti(:,1);
            tv(yindnn) = ti(:,2);
            sv(xindnn) = si(:,1);
            sv(yindnn) = si(:,2);
            last = next + 1;
        end
        
        % normalize tv and sv across all vertices
        tvnorm = sqrt(sum(tv.^2));
        svnorm = sqrt(sum(sv.^2));
        
        tv = tv./tvnorm;
        sv = sv./svnorm;
        
        % compute dot product of each constructed direction with
        % eigenvector
        t_proj = zeros(NDOF,1);
        s_proj = zeros(NDOF,1);
        for dd = 1:NDOF
            t_proj(dd) = sum(abs(tv.*evectmp(:,dd)));
            s_proj(dd) = sum(abs(sv.*evectmp(:,dd)));
        end
        
        % save projection
        projections{pp,4} = t_proj;
        projections{pp,5} = s_proj;
    end


    fprintf('\t ** Finishing storing simulation data...\n');

    % number of quartic modes
    lambdaMatList{ss}       = lambdaMatrix;
    stiffEVMatList{ss}      = stiffEVMatrix;
    NqList{ss}              = nq;
    vertexPrList{ss}        = vertexPR;
    cellPrList{ss}          = cellPR;
    projList{ss}            = projections;
    
    
    
    % -- cell energy data
    jFrameList{ss}          = jframes;
    
    enStress                = cellEnergyData.stress;
    enPhi                   = cellEnergyData.phi;
    enP                     = 0.5*(enStress(:,1) + enStress(:,3))/(NCELLS*Lx*Ly);
    enCalA                  = cellEnergyData.calA;
    enCalA0                 = cellEnergyData.calA0;
    enPList{ss}             = enP;
    enPhiList{ss}           = enPhi;
    enCalAList{ss,1}        = enCalA;
    enCalAList{ss,2}        = enCalA0;
end

% delete any entries that were skipped
NCELLSList(simSkip)         = [];
NDOFList(simSkip)           = [];
NvList(simSkip)             = [];
LList(simSkip,:)            = [];
a0List(simSkip)             = [];
l0List(simSkip)             = [];
pList(simSkip)              = [];
NvvList(simSkip)            = [];
NccList(simSkip)            = [];
NqList(simSkip)             = [];
lambdaMatList(simSkip)      = [];
stiffEVMatList(simSkip)     = [];
vertexPrList(simSkip)       = [];
cellPrList(simSkip)         = [];
projList(simSkip)           = [];
jFrameList(simSkip)         = [];
enPList(simSkip)            = [];
enPhiList(simSkip)          = [];
enCalAList(simSkip,:)         = [];


% total number of sims after removing empties
NSIMS = sum(~simSkip);


%% Save to matfile

fprintf('\t ** saving to %s\n',saveStr);

save(saveStr,...
    'NSIMS','simStr','NCELLSList','NDOFList','NvList','LList','a0List','l0List',...
    'NFRAMEList','pList','NvvList','NccList','NqList',...
    'lambdaMatList','stiffEVMatList','vertexPrList','cellPrList','projList',...
    'jFrameList','enPList','enPhiList','enCalAList');

end