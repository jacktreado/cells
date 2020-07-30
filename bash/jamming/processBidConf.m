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
riList          = cell(NSIM,1);         % number of rattlers list

% containers for extracted VDOS data
pList           = cell(NSIM,1);         % list of positive pressures
NvvList         = cell(NSIM,1);         % # of vertex-vertex contacts
NccList         = cell(NSIM,1);         % # of cell-cell contacts
NvvrList        = cell(NSIM,1);         % # of vertex-vertex contacts after rattler removal
NccrList        = cell(NSIM,1);         % # of cell-cell contacts after rattler removal
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
    
    % get mechanical parameters
    finfo = dir(jamfstr);
    fname = finfo.name;
    calA0Str = sscanf(fname,['bidconf_N' num2str(NCELLS) '_NV' num2str(nv(1)) '_calA%[0-9.]_kl']);
    Klstr = sscanf(fname,['bidconf_N' num2str(NCELLS) '_NV' num2str(nv(1)) '_calA' calA0Str '_kl%[0-9.]_']);
    Kbstr = sscanf(fname,['bidconf_N' num2str(NCELLS) '_NV' num2str(nv(1)) '_calA' calA0Str '_kl' Klstr '_kb%[0-9.]_']);
    
    Kl = str2double(Klstr);
    Kb = str2double(Kbstr);
    
    % compute matrix to map cell-and-vertex indices to global indices
    maxNV   = max(nv);
    vvi     = zeros(NCELLS,maxNV);
    jj      = 1;
    for ii = 1:NCELLS
        for vv = 1:nv(ii)
            vvi(ii,vv) = jj;
            jj = jj + 1;
        end
    end
        
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
            evtmp = evectmp(:,dd);
            last = 0;
            for nn = 1:NCELLS
                for vv = 1:nv(nn)
                    xi = 2*vv - 1 + last;
                    yi = 2*vv + last;
                    t_proj(dd) = t_proj(dd) + abs(tv(xi)*evtmp(xi) + tv(yi)*evtmp(yi));
                    s_proj(dd) = s_proj(dd) + abs(sv(xi)*evtmp(xi) + sv(yi)*evtmp(yi));
                end
                last = last + 2*nv(nn);
            end
        end
        
        % save projection
        projections{pp,4} = t_proj;
        projections{pp,5} = s_proj;
    end

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
    
    
    
    % -- compute number of rattlers
    fprintf('\t ** -- Computing rattler data...\n');
    
    cxpos = cellJamData.cxpos;
    cypos = cellJamData.cypos;
    
    ritmp   = zeros(NPOSFRAMES,NCELLS);
    spdof   = zeros(NPOSFRAMES,NCELLS);
    Nvvr    = zeros(NPOSFRAMES,1);
    Nccr    = zeros(NPOSFRAMES,1);
    for pp = 1:NPOSFRAMES
        nvvMeas = 0;
        nccMeas = 0;
        
        % contact matrices
        cij = zeros(NCELLS);
        bcij = zeros(NCELLS);
        
        % vertex positions and sizes
        xtmp = xpos(pp,:);
        ytmp = ypos(pp,:);
        l0tmp = l0Pos(pp,:);
        
        xall = zeros(NVTOT,1);
        yall = zeros(NVTOT,1);
        
        % contacts between each vertex
        cab = false(NVTOT);
        
        % loop over vertex pairs, try to get correctly Nvv and Ncc
        for ii = 1:NCELLS
            xii = xtmp{ii};
            yii = ytmp{ii};
            dii = l0tmp(ii);
            
            xall(vvi(ii,1):vvi(ii,nv(ii))) = xii;
            yall(vvi(ii,1):vvi(ii,nv(ii))) = yii;
            
            for jj = ii+1:NCELLS
                xjj = xtmp{jj};
                yjj = ytmp{jj};
                djj = l0tmp(jj);
                
                ccFound = 0;
                
                for aa = 1:nv(ii)
                    for bb = 1:nv(jj)
                        sab = 0.5*(dii + djj);
                        
                        % scale contact distance by a little bit to make
                        % contacts easier to count
                        sab = sab*(1.0 + 1e-8);
                        
                        dx = xjj(bb) - xii(aa);
                        dx = dx - Lx*round(dx/Lx);
                        
                        dy = yjj(bb) - yii(aa);
                        dy = dy - Ly*round(dy/Ly);
                        
                        dr = sqrt(dx*dx + dy*dy);
                        
                        % if vertices overlap, increase contact count
                        if dr < sab
                            % add contacts
                            if ccFound == 0
                                nccMeas = nccMeas + 1;
                                bcij(ii,jj) = 1;
                                bcij(jj,ii) = 1;
                                ccFound = 1;
                            end
                            nvvMeas = nvvMeas + 1;
                            cij(ii,jj) = cij(ii,jj) + 1;
                            cij(jj,ii) = cij(jj,ii) + 1;
                            
                            % save indices of contacting vertices
                            aInd = vvi(ii,aa);
                            bInd = vvi(jj,bb);
                            cab(aInd,bInd) = true;
                            cab(bInd,aInd) = true;
                        end
                    end
                end
            end
        end
        
        dnvv = NvvPos(pp) - nvvMeas;
        dncc = NccPos(pp) - nccMeas;
        fprintf('\t\t ** pp = %d, dnvv = %d, dncc = %d\n',pp,dnvv,dncc);
        
        % get initial cell contacts
        zcc = sum(bcij,1);
        
        % number of rattlers, marginal rattlers
        nm = 1;
        it = 0;
        itmax = 1e3;
        while (nm > 0 && it < itmax)
            
            % reset marginal count to 0
            nm = 0;
            
            % remove all cells with only 1 cell-cell contact
            ct1_inds = (zcc == 1);
            nct = sum(ct1_inds);
            
            nm = nm + nct;
            
            cij(ct1_inds,:) = zeros(nct,NCELLS);
            cij(:,ct1_inds) = zeros(NCELLS,nct);
            
            bcij(ct1_inds,:) = zeros(nct,NCELLS);
            bcij(:,ct1_inds) = zeros(NCELLS,nct);
            
            % update zcc
            zcc = sum(bcij,1);
            
            % remove all cells with 2 cell-cell contacts, at least 1 is
            % marginal
            mvvcts = sum(cij == 1,2);
            ct2_inds = ((mvvcts == 1 | mvvcts == 2) & (zcc == 2)');
            nct = sum(ct2_inds);
            
            nm = nm + nct;
            
            cij(ct2_inds,:) = zeros(nct,NCELLS);
            cij(:,ct2_inds) = zeros(NCELLS,nct);
            
            bcij(ct2_inds,:) = zeros(nct,NCELLS);
            bcij(:,ct2_inds) = zeros(NCELLS,nct);
            
            % update zcc
            zcc = sum(bcij,1);
            
            % increase iterate check
            it = it + 1;
        end
        if it == itmax
            error('Error in processBidConf, iteration max found while removing rattlers.');
        end
        
        % tally up rattlers
        ri = (zcc == 0);
        
        
        % compute nvv and ncc after rattler removal
        zvv = sum(cij,1);
        
        Nvvr(pp) = NCELLS*mean(zvv)/2;
        Nccr(pp) = NCELLS*mean(zcc)/2;
        
        % save for this pressure
        ritmp(pp,:) = ri';
        
        % -- compute VDOS for single particles, save effective d.o.f.
        cxtmp = cxpos(pp,:);
        cytmp = cypos(pp,:);
        l0tmp = l0Pos(pp,:);
        a0tmp = a0Pos(pp,:);
        for ii = 1:NCELLS
            % vertex coordinates for cell ii
            vx = xtmp{ii};
            vy = ytmp{ii};
            
            % get relative coordinates
            rx = vx - cxtmp(ii);
            rx = rx - Lx*round(rx./Lx);
            
            ry = vy - cytmp(ii);
            ry = ry - Ly*round(ry./Ly);

            % get pinned vertex positions and sizes
            px = zeros(nv(ii)*100,1);
            py = zeros(nv(ii)*100,1);
            l0p = zeros(nv(ii)*100,1);
            last = 1;
            for vv = 1:nv(ii)
                vInd = vvi(ii,vv);
                cInds = cab(:,vInd);
                NVC = sum(cInds);
                
                next = last + NVC - 1;
                px(last:next) = xall(cInds);
                py(last:next) = yall(cInds);
                
                % get sizes
                smallInds = find(cInds) < 0.5*NVTOT;
                possInds = 1:NVC;
                smallLocs = last + possInds(smallInds) - 1;
                largeLocs = last + possInds(~smallInds) - 1;
                l0p(smallLocs) = repmat(l0tmp(1),length(smallLocs),1);
                l0p(largeLocs) = repmat(l0tmp(end),length(largeLocs),1);
                
                % increment
                last = next + 1;
            end
            px(last:end) = [];
            py(last:end) = [];
            l0p(last:end) = [];
            
            [pxUQ, IA, ~] = unique(px);
            pyUQ = py(IA);
            
            % call VDOS function on this cell with contacts
            calA0ii = (nv(ii)*l0tmp(ii))^2/(4.0*pi*a0tmp(ii));
            [Ha, Sa, Hl, Sl, Hb, Sb, Hvv, Svv] = dpVDOSExternalForces(vx, vy, rx, ry, Lx, Ly, pxUQ, pyUQ, a0tmp(ii), l0p, Kl, Kb, calA0ii);
            
            % get matrices
            H = Ha + Hl + Hb + Hvv;
            S = Sa + Sl + Sb + Svv;
            D = H + S;
            
            if zvv(ii) == 0
                
                
                [singleEV, devals] = eig(D);
                lambda = diag(devals);
                

                figure(1), clf, hold on, box on;
                for jj = 1:NCELLS
                    xv = xtmp{jj};
                    yv = ytmp{jj};
                    ww = l0tmp(jj);  
                    for vv = 1:nv(jj)
                        xx = xv(vv);
                        yy = yv(vv);
                        if ri(jj) == 0
                            rectangle('Position',[xx-0.5*ww,yy-0.5*ww,ww,ww],'Curvature',[1 1],'EdgeColor','k');
                        else
                            rectangle('Position',[xx-0.5*ww,yy-0.5*ww,ww,ww],'Curvature',[1 1],'EdgeColor','r');
                        end
                    end
                end
                quiver(vx,vy,singleEV(1:2:(2*nv(ii)-1),1),singleEV(2:2:(2*nv(ii)),1),'-r','linewidth',2);
                plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-','linewidth',1.5);
                axis equal;
                ax = gca;
                ax.XTick = [];
                ax.YTick = [];
                ax.XLim = [-0.25 1.25].*Lx;
                ax.YLim = [-0.25 1.25].*Ly;
                
                test = 1;
            end
        end
        
    end
    
    % save rattker matrix to list
    riList{ss} = ritmp;
    NvvrList{ss} = Nvvr;
    NccrList{ss} = Nccr;
    
    fprintf('\t ** Finishing storing simulation data for ss = %d...\n',ss);
end

% delete any entries that were skipped
simList(simSkip)            = [];
NCELLSList(simSkip)         = [];
NDOFList(simSkip)           = [];
NvList(simSkip)             = [];
LList(simSkip,:)            = [];
a0List(simSkip)             = [];
l0List(simSkip)             = [];
pList(simSkip)              = [];
NvvList(simSkip)            = [];
NccList(simSkip)            = [];
NvvrList(simSkip)           = [];
NccrList(simSkip)           = [];
NqList(simSkip)             = [];
lambdaMatList(simSkip)      = [];
stiffEVMatList(simSkip)     = [];
vertexPrList(simSkip)       = [];
cellPrList(simSkip)         = [];
projList(simSkip)           = [];
jFrameList(simSkip)         = [];
enPList(simSkip)            = [];
enPhiList(simSkip)          = [];
enCalAList(simSkip,:)       = [];
riList(simSkip)             = [];

% total number of sims after removing empties
NSIMS = sum(~simSkip);


%% Save to matfile

fprintf('\t ** saving to %s\n',saveStr);

save(saveStr,...
    'NSIMS','simList','NCELLSList','NDOFList','NvList','LList','a0List','l0List',...
    'NFRAMEList','pList','NvvList','NccList','NvvrList','NccrList','NqList',...
    'lambdaMatList','stiffEVMatList','vertexPrList','cellPrList','projList',...
    'jFrameList','enPList','enPhiList','enCalAList','riList');

end












%% FUNCTION to compute VDOS for single particle under external forces

function [Ha, Sa, Hl, Sl, Hb, Sb, Hvv, Svv] = dpVDOSExternalForces(vx, vy, rx, ry, Lx, Ly, px, py, a0, l0p, Kl, Kb, calA0)

% get number of vertices (check inputs)
NVx = size(rx,1);
NVy = size(ry,1);
if (NVx ~= NVy)
    fprintf('lengths of vx and vy do not match in dpVDOS function...\n');
    fprintf('size(vx) = %d %d, size(vy) = %d %d, ending\n',size(rx,1),size(rx,2),size(ry,1), size(ry,2));
    error('size of inputs is incorrect');
else
    NV = NVx;
end

% get polygon area and perimeter
a = area(vx, vy);

% get polygon perimeter and segment lengths
[~, l] = perimeter(vx, vy);

% l0 bar
l0 = sqrt(4*pi*calA0*a0)/NV;

%% Calculate dynamical matrix (stiffness and stress matrices independently)

% stress matrix magnitudes
delA = (a/a0) - 1;
dell = zeros(NV,1);
for vv = 1:NV
    dell(vv) = (l(vv)/l0) - 1;
end

% bending energy constant
eb = (Kb*NV*calA0)/(4*pi^2);
fb = eb/l0^2;

% perimeter spring constant 
kl = Kl/l0;

% number of degrees of freedom
dof = 2*NV;

% area term matrices
Ha = zeros(dof);
Sa = zeros(dof);

Hl = zeros(dof);
Sl = zeros(dof);

Hb = zeros(dof);
Sb = zeros(dof);


% loop over vertices, calculate each DM matrix element
for r = 1:NV
    
    % wrap around ring of vertices
    rp2 = r + 2;
    rp1 = r + 1;
    rm1 = r - 1;
    rm2 = r - 2;
    if r == 2
        rm2 = NV;
    elseif r == 1
        rm1 = NV;
        rm2 = NV-1;
    elseif r == NV
        rp1 = 1;
        rp2 = 2;
    elseif r == NV-1
        rp2 = 1;
    end
    
    % dof elements
    kx          = 2*r - 1;
    ky          = 2*r;
    
    kxp1        = 2*rp1 - 1;
    kyp1        = 2*rp1;
    
    kxp2        = 2*rp2 - 1;
    kyp2        = 2*rp2;
    
    % store coordinates
    xrp2        = rx(rp2);
    xrp1        = rx(rp1);
    xr          = rx(r);
    xrm1        = rx(rm1);
    xrm2        = rx(rm2);
    
    yrp2        = ry(rp2);
    yrp1        = ry(rp1);
    yr          = ry(r);
    yrm1        = ry(rm1);
    yrm2        = ry(rm2);
    
    % segment length vectors
    lrm2x       = xrm1 - xrm2;
    lrm2y       = yrm1 - yrm2;
    
    lrm1x       = xr - xrm1;
    lrm1y       = yr - yrm1;
    
    lrx         = xrp1 - xr;
    lry         = yrp1 - yr;
    
    lrp1x       = xrp2 - xrp1;
    lrp1y       = yrp2 - yrp1;
    
    % segment lengths
    lrm1        = l(rm1);
    lr          = l(r);
    
    % PERIMETER ENERGY: compute derivatives of segment lengths w.r.t. coordinates
    
    % derivatives of lrm1
    dlrm1_dxr   = lrm1x/lrm1;
    dlrm1_dyr   = lrm1y/lrm1;
    
    % derivatives of lr
    dlr_dxrp1   = lrx/lr;
    dlr_dyrp1   = lry/lr;
    
    dlr_dxr     = -dlr_dxrp1;
    dlr_dyr     = -dlr_dyrp1;
    
    
    % -- stiffness matrix
    
    % main diagonal
    Hl(kx,kx)       = kl*(dlrm1_dxr*dlrm1_dxr + dlr_dxr*dlr_dxr);
    Hl(ky,ky)       = kl*(dlrm1_dyr*dlrm1_dyr + dlr_dyr*dlr_dyr);
    
    Hl(kx,ky)       = kl*(dlrm1_dxr*dlrm1_dyr + dlr_dxr*dlr_dyr);
    Hl(ky,kx)       = Hl(kx,ky);
    
    % 1off diagonal
    Hl(kx,kxp1)     = kl*dlr_dxr*dlr_dxrp1;
    Hl(ky,kyp1)     = kl*dlr_dyr*dlr_dyrp1;
    
    Hl(kx,kyp1)     = kl*dlr_dxr*dlr_dyrp1;
    Hl(ky,kxp1)     = kl*dlr_dyr*dlr_dxrp1;
    
    % enforce symmetry in lower triangle
    Hl(kxp1,kx)     = Hl(kx,kxp1);
    Hl(kyp1,ky)     = Hl(ky,kyp1);
    
    Hl(kyp1,kx)     = Hl(kx,kyp1);
    Hl(kxp1,ky)     = Hl(ky,kxp1);
    
    
    % -- stress matrix elements
    
    % main diagonal
    Sl(kx,kx)       = Kl*((dell(rm1)/lrm1)*(1 - (dlrm1_dxr*dlrm1_dxr)) + (dell(r)/lr)*(1 - (dlr_dxr*dlr_dxr)));
    Sl(ky,ky)       = Kl*((dell(rm1)/lrm1)*(1 - (dlrm1_dyr*dlrm1_dyr)) + (dell(r)/lr)*(1 - (dlr_dyr*dlr_dyr)));
    
    Sl(kx,ky)       = -Kl*((dell(rm1)/lrm1)*dlrm1_dxr*dlrm1_dyr + (dell(r)/lr)*dlr_dxr*dlr_dyr);
    Sl(ky,kx)       = Sl(kx,ky);
    
    % 1off diagonal
    Sl(kx,kxp1)     = Kl*(dell(r)/lr)*((dlr_dxrp1*dlr_dxrp1) - 1);
    Sl(ky,kyp1)     = Kl*(dell(r)/lr)*((dlr_dyrp1*dlr_dyrp1) - 1);
    
    Sl(kx,kyp1)     = Kl*(dell(r)/lr)*dlr_dxrp1*dlr_dyrp1;
    Sl(ky,kxp1)     = Kl*(dell(r)/lr)*dlr_dyrp1*dlr_dxrp1;
    
    % enforce symmetry in lower triangle
    Sl(kxp1,kx)     = Sl(kx,kxp1);
    Sl(kyp1,ky)     = Sl(ky,kyp1);
    
    Sl(kyp1,kx)     = Sl(kx,kyp1);
    Sl(kxp1,ky)     = Sl(ky,kxp1);
    
    
    
    % BENDING ENERGY: compute derivatives of curvature w.r.t. coordinates
    
    
    % -- dimensionless curvatures
    kaprm1          = sqrt((lrm1x - lrm2x)^2 + (lrm1y - lrm2y)^2)/l0;
    kapr            = sqrt((lrx - lrm1x)^2 + (lry - lrm1y)^2)/l0;
    kaprp1          = sqrt((lrp1x - lrx)^2 + (lrp1y - lry)^2)/l0;

    % -- curvature derivatives
    
    % derivatives of kaprm1
    dkaprm1_dxr     = (lrm1x - lrm2x)/(kaprm1*l0^2);
    dkaprm1_dyr     = (lrm1y - lrm2y)/(kaprm1*l0^2);
    
    % derivatives of kapr
    dkapr_dxrp1     = (lrx - lrm1x)/(kapr*l0^2);
    dkapr_dyrp1     = (lry - lrm1y)/(kapr*l0^2);
    dkapr_dxr       = -2*dkapr_dxrp1;
    dkapr_dyr       = -2*dkapr_dyrp1;
    
    % derivatives of kaprp1
    dkaprp1_dxr     = (lrp1x - lrx)/(kaprp1*l0^2);
    dkaprp1_dyr     = (lrp1y - lry)/(kaprp1*l0^2);
    dkaprp1_dxrp1   = -2*dkaprp1_dxr;
    dkaprp1_dyrp1   = -2*dkaprp1_dyr;
    dkaprp1_dxrp2   = dkaprp1_dxr;
    dkaprp1_dyrp2   = dkaprp1_dyr;
    
    
    
    
    % BENDING ENERGY: stiffness and stress matrices
    
    
    % -- stiffness matrix elements
    
    % block-diagonal terms
    Hb(kx,kx)       = eb*(dkaprm1_dxr^2 + dkapr_dxr^2 + dkaprp1_dxr^2);
    Hb(ky,ky)       = eb*(dkaprm1_dyr^2 + dkapr_dyr^2 + dkaprp1_dyr^2);
    
    Hb(kx,ky)       = eb*(dkaprm1_dxr*dkaprm1_dyr + dkapr_dxr*dkapr_dyr + dkaprp1_dxr*dkaprp1_dyr);
    Hb(ky,kx)       = Hb(kx,ky);
    
    % 1off block-diagonal terms
    Hb(kx,kxp1)     = eb*(dkapr_dxr*dkapr_dxrp1 + dkaprp1_dxr*dkaprp1_dxrp1);
    Hb(ky,kyp1)     = eb*(dkapr_dyr*dkapr_dyrp1 + dkaprp1_dyr*dkaprp1_dyrp1);
    
    Hb(kx,kyp1)     = eb*(dkapr_dxr*dkapr_dyrp1 + dkaprp1_dxr*dkaprp1_dyrp1);
    Hb(ky,kxp1)     = eb*(dkapr_dyr*dkapr_dxrp1 + dkaprp1_dyr*dkaprp1_dxrp1);
    
    % 2off block-diagonal terms
    Hb(kx,kxp2)     = eb*dkaprp1_dxr*dkaprp1_dxrp2;
    Hb(ky,kyp2)     = eb*dkaprp1_dyr*dkaprp1_dyrp2;
    
    Hb(kx,kyp2)     = eb*dkaprp1_dxr*dkaprp1_dyrp2;
    Hb(ky,kxp2)     = eb*dkaprp1_dyr*dkaprp1_dxrp2;
    
    % enforce symmetry in lower triangle
    Hb(kxp1,kx)     = Hb(kx,kxp1);
    Hb(kyp1,ky)     = Hb(ky,kyp1);
    
    Hb(kxp1,ky)     = Hb(ky,kxp1);
    Hb(kyp1,kx)     = Hb(kx,kyp1);
    
    Hb(kxp2,kx)     = Hb(kx,kxp2);
    Hb(kyp2,ky)     = Hb(ky,kyp2);
    
    Hb(kyp2,kx)     = Hb(kx,kyp2);
    Hb(kxp2,ky)     = Hb(ky,kxp2);
    
    
    % -- stress matrix elements
    
    % block diagonal
    Sb(kx,kx)       = fb*(6 - (l0*dkaprm1_dxr)^2 - (l0*dkapr_dxr)^2 - (l0*dkaprp1_dxr)^2);
    Sb(ky,ky)       = fb*(6 - (l0*dkaprm1_dyr)^2 - (l0*dkapr_dyr)^2 - (l0*dkaprp1_dyr)^2);
    
    Sb(kx,ky)       = -eb*(dkaprm1_dxr*dkaprm1_dyr + dkapr_dxr*dkapr_dyr + dkaprp1_dxr*dkaprp1_dyr);
    Sb(ky,kx)       = Sb(kx,ky);
    
    % 1off block diagonal
    Sb(kx,kxp1)     = -2*fb*(2 - (l0*dkapr_dxrp1)^2 - (l0*dkaprp1_dxr)^2);
    Sb(ky,kyp1)     = -2*fb*(2 - (l0*dkapr_dyrp1)^2 - (l0*dkaprp1_dyr)^2);
    
    Sb(kx,kyp1)     = -eb*(dkapr_dxr*dkapr_dyrp1 + dkaprp1_dxr*dkaprp1_dyrp1);
    Sb(ky,kxp1)     = -eb*(dkapr_dyr*dkapr_dxrp1 + dkaprp1_dyr*dkaprp1_dxrp1);
    
    % 2off block diagonal
    Sb(kx,kxp2)     = fb*(1 - (l0*dkaprp1_dxr)^2);
    Sb(ky,kyp2)     = fb*(1 - (l0*dkaprp1_dyr)^2);
    
    Sb(kx,kyp2)     = -eb*dkaprp1_dxr*dkaprp1_dyrp2;
    Sb(ky,kxp2)     = -eb*dkaprp1_dyr*dkaprp1_dxrp2;
    
    % enforce symmetry in lower triangle
    Sb(kxp1,kx)     = Sb(kx,kxp1);
    Sb(kyp1,ky)     = Sb(ky,kyp1);
    
    Sb(kxp1,ky)     = Sb(ky,kxp1);
    Sb(kyp1,kx)     = Sb(kx,kyp1);
    
    Sb(kxp2,kx)     = Sb(kx,kxp2);
    Sb(kyp2,ky)     = Sb(ky,kyp2);
    
    Sb(kxp2,ky)     = Sb(ky,kxp2);
    Sb(kyp2,kx)     = Sb(kx,kyp2);
    
    % AREA ELASTICITY: compute derivatives of area (a) w.r.t. xr and yr
    da_dxr      = 0.5*(yrm1 - yrp1);
    da_dyr      = 0.5*(xrp1 - xrm1);
    
    % area elasticity stress matrix elements
    Sa(kx,kyp1) = 0.5*delA;
    Sa(ky,kxp1) = -0.5*delA;

    Sa(kyp1,kx) = Sa(kx,kyp1);
    Sa(kxp1,ky) = Sa(ky,kxp1);
    
    % loop over vertices s = r to NV (keep r to get block diagonal
    % elements) just for area elasticity
    for s = r:NV
        
        % wrap sp1 and sm1 around ring of vertices
        sp1 = s + 1;
        sm1 = s - 1;
        if s == 1
            sm1 = NV;
        elseif s == NV
            sp1 = 1;
        end
        
        % dof elements
        lx       = 2*s - 1;
        ly       = 2*s;

        % store coordinates
        xsp1        = rx(sp1);
        xsm1        = rx(sm1);

        ysp1        = ry(sp1);
        ysm1        = ry(sm1);

        % AREA ELASTICITY: compute derivatives of area (a) w.r.t. xr and yr
        da_dxs      = 0.5*(ysm1 - ysp1);
        da_dys      = 0.5*(xsp1 - xsm1);
        
        % COMPUTE MATRIX ELEMENTS
        
        % area elasticity stiffness matrix elements 
        Ha(kx,lx) = da_dxr*da_dxs;
        Ha(kx,ly) = da_dxr*da_dys;
        
        Ha(ky,lx) = da_dyr*da_dxs;
        Ha(ky,ly) = da_dyr*da_dys;
        
        Ha(lx,kx) = Ha(kx,lx);
        Ha(ly,kx) = Ha(kx,ly);
        
        Ha(lx,ky) = Ha(ky,lx);
        Ha(ly,ky) = Ha(ky,ly);
    end
end


% compute dynamical matrix element based on location of external pinned
% vertices
Hvv = zeros(dof);
Svv = zeros(dof);
NPINS = length(px);
for pp = 1:NPINS
    % get coordinates of external vertices
    xtmp = px(pp);
    ytmp = py(pp);

    % get size of external vertices
    sp = l0p(pp);
    
    % loop over internal vertices, compute forces
    for vv = 1:NV
        % get coordinates of internal vertex
        xv = vx(vv);
        yv = vy(vv);
        
        % get contact size
        sij = 0.5*(l0 + sp);

        % get distances to closest minimum
        dx = xtmp - xv;
        dx = dx - Lx*round(dx/Lx);
        
        dy = ytmp - yv;
        dy = dy - Ly*round(dy/Ly);
        
        dr = sqrt(dx*dx + dy*dy);
        
        % matrix elements
        xi = 2*vv - 1;
        yi = 2*vv;

        if dr < l0
            % spring constant
            kij = 1/(sij*dr);
            
            % overlap
            h = dr/sij;
            
            % derivatives
            dr_dxi = -dx/dr;
            dr_dyi = -dy/dr;
            
            % stiffness matrix elements
            Hvv(xi,xi) = Hvv(xi,xi) + dr_dxi*dr_dxi;
            Hvv(yi,yi) = Hvv(yi,yi) + dr_dyi*dr_dyi;
            
            Hvv(xi,yi) = Hvv(xi,yi) + dr_dxi*dr_dyi;
            Hvv(yi,xi) = Hvv(yi,xi) + dr_dyi*dr_dxi;
            
            % stress matrix elements
            Svv(xi,xi) = Svv(xi,xi) - kij*(1.0 - h)*(dr_dyi*dr_dyi);
            Svv(yi,yi) = Svv(xi,xi) - kij*(1.0 - h)*(dr_dxi*dr_dxi);
            
            Svv(xi,yi) = Svv(xi,yi) + kij*(1.0 - h)*(dr_dxi*dr_dyi);
            Svv(yi,xi) = Svv(yi,xi) + kij*(1.0 - h)*(dr_dxi*dr_dyi);
        end
    end
end


end




function [P, lens] = perimeter(vx, vy)
%% FUNCTION to calculate area of a polygonal DPM particle based of vertex coordinates

% check inputs
if sum(size(vx)) ~= sum(size(vy))
    fprintf('In area function, size vx = %d %d, size vy = %d %d\n',size(vx,1),size(vx,2),size(vy,1),size(vy,2));
    error('Sizes of inputs vx and vy do not match, ending');
end

% get number of vertices
NV = length(vx);

% loop indices
loopi = [1:NV 1];

% perimeter sum
dx      = vx(loopi(2:end)) - vx(loopi(1:end-1));
dy      = vy(loopi(2:end)) - vy(loopi(1:end-1));
lens    = sqrt(dx.^2 + dy.^2);
P       = sum(lens);

end


function A = area(vx,vy)
%% FUNCTION to calculate area of a polygonal DPM particle based of vertex coordinates

% check inputs
if sum(size(vx)) ~= sum(size(vy))
    fprintf('In area function, size vx = %d %d, size vy = %d %d\n',size(vx,1),size(vx,2),size(vy,1),size(vy,2));
    error('Sizes of inputs vx and vy do not match, ending');
end

% get number of vertices
NV = length(vx);

%% Calculate area

% get vertex with max y values
[~, kStart] = max(vy);

% get left and right of kStart
km1  = mod((kStart-1) + NV - 1, NV) + 1;
kp1  = mod((kStart-1) + 1, NV) + 1;

% find out what direction is to the right of kStart
if vx(kp1) > vx(km1)
    cwdir = true;
elseif vx(kp1) < vx(km1)
    cwdir = false;
else
    km2 = mod((km1-1) + NV - 1, NV) + 1;
    kp2 = mod((kp1-1) + 1,NV) + 1;
    if vx(kp2) >= vx(km2)
        cwdir = true;
    else
        cwdir = false;
    end
end

% set up indices
kcurr = kStart;

% initialize area
polyArea = 0.0;

% loop over vertices
cnt = 0;
while cnt < NV
    % increment counter
    cnt = cnt + 1;
    
    % determine index of next vertex
    if cwdir
        knext = mod((kcurr-1) + 1,NV) + 1;
    else
        knext = mod((kcurr-1) + NV - 1, NV) + 1;
    end
    
    % get coordinates
    xcurr = vx(kcurr);
    ycurr = vy(kcurr);
    
    xnext = vx(knext);
    ynext = vy(knext);
    
    % check Dong's version
    polyArea = polyArea + 0.5*(xnext*ycurr - xcurr*ynext);
    
    % update kcurr for next iteration
    kcurr = knext;
end

% halve area for return
A = polyArea;

end

