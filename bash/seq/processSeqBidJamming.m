function processSeqBidJamming(simStr, saveStr)
%% FUNCTION to read in vdos files and jamming file, process and save data to matfiles

% get list of files
flist = dir([simStr '*.vdos']);
NSIM = length(flist);
if NSIM == 0
    fprintf('No files found using str %s\n',simStr);
    error('No files found, ending.');
end

%% Loop over files, save data

% simulation arrays
simSkip             = false(NSIM,1);        % whether or not to skip file in list
simList             = cell(NSIM,1);         % simulation file name list
NCELLSList          = zeros(NSIM,1);        % number of cells in each sim
NvList              = cell(NSIM,1);         % # of vertices on each particle
LList               = zeros(NSIM,2);        % box lengths

% particle arrays
zcList              = cell(NSIM,1);         % cc contacts per particle
zvList              = cell(NSIM,1);         % vv contacts per particle
a0List              = cell(NSIM,1);         % particle preferred area list
l0List              = cell(NSIM,1);         % vertex size list
calAList            = cell(NSIM,2);         % instaneous shape parameter of each particle at jamming
aList               = cell(NSIM,1);
phiJList            = cell(NSIM,2);         % packing fractions (both from a and a0) at jamming

% voronoi arrays
voroAreasList       = cell(NSIM,1);         % voronoi areas
voroCalAList        = cell(NSIM,1);         % voronoi calA

% VDOS info
evalsList           = cell(NSIM,1);         % list of eigenvalues at jamming
hvalsList           = cell(NSIM,1);         % list of projections of eigenvectors onto H
pratioList          = cell(NSIM,1);         % list of participation ratios

%% Loop over simulations


for ss = 1:NSIM
    % file information
    dirname                 = flist(ss).folder;
    vdosfname               = flist(ss).name;
    jamfname                = [vdosfname(1:end-5) '.jam'];
    
    jamfstr                 = [dirname '/' jamfname];
    vdosfstr                = [dirname '/' vdosfname];
    
    % test to make sure that all files exist and have data
    if exist(jamfstr,'file') > 0
        jamfinfo = dir(jamfstr);
        if jamfinfo.bytes == 0
            fprintf('jam file %s exists but is empty, skipping sim...\n',jamfstr);
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
                    simList{ss} = vdosfname;
                end
            else
                fprintf('vdos file %s does not exist in location, skipping sim...\n',vdosfstr);
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
    
    % get position data for quenched packing
    quenchedPackingData = readQuenchedPackings(jamfstr);
    
    % parse position data
    NCELLS      = quenchedPackingData.NCELLS;
    L           = quenchedPackingData.L;
    nv          = quenchedPackingData.nv(1,:);
    xpos        = quenchedPackingData.xpos(1,:);
    ypos        = quenchedPackingData.ypos(1,:);
    zc          = quenchedPackingData.zc(1,:);
    zv          = quenchedPackingData.zv(1,:);
    l0          = quenchedPackingData.l0(1,:);
    a0          = quenchedPackingData.a0(1,:);
    
    % compute perimeter, area, shape parameters and packing fractions
    a = zeros(NCELLS,1);
    p = zeros(NCELLS,1);
    calA = zeros(NCELLS,1);
    calA0 = zeros(NCELLS,1);
    phiJ = 0.0;
    phi0J = 0.0;
    for nn = 1:NCELLS
        x = xpos{nn};
        y = ypos{nn};
        
        dx = x([2:end 1]) - x;
        dy = y([2:end 1]) - y;
        l = sqrt(dx.^2 + dy.^2);
        
        a(nn) = polyarea(x,y);
        p(nn) = sum(l);
        calA(nn) = p(nn)^2/(4.0*pi*a(nn));
        calA0(nn) = (nv(nn)*l0(nn))^2/(4.0*pi*a0(nn));
        
        phiJ = phiJ + a(nn) + 0.25*pi*l0(nn)^2*(0.5*nv(nn) - 1);
        phi0J = phi0J + a0(nn) + 0.25*pi*l0(nn)^2*(0.5*nv(nn) - 1);
    end
    phiJ = phiJ/(L(1)*L(2));
    phi0J = phi0J/(L(1)*L(2));
    
    % save particle information to lists
    NCELLSList(ss)      = NCELLS;
    NvList{ss}          = nv;
    LList(ss,1)         = L(1);
    LList(ss,2)         = L(2);
    
    zcList{ss}          = zc;
    zvList{ss}          = zv;
    a0List{ss}          = a0;
    l0List{ss}          = l0;
    
    calAList{ss,1}      = calA;
    calAList{ss,2}      = calA0;
    
    phiJList{ss,1}      = phiJ;
    phiJList{ss,2}      = phi0J;
    
    
    % read in vdos data
    fid = fopen(vdosfstr);

    % dof
    dof = textscan(fid,'%f',1);
    dof = dof{1};

    % eigenvalue data
    evals = textscan(fid,'%f',dof);
    evals = evals{1};


    % stiffness eval data
    hvals = textscan(fid,'%f',dof);
    hvals = hvals{1};
    
    % check end of file, it not finished then
    % assume rest of data is vertex participation
    % ratio
    % check for participation ratio
    if ~feof(fid)
        prtmp = textscan(fid,'%f',dof);
        pr = prtmp{1};
        pratioList{ss} = pr;
    end

    % close the file
    fclose(fid);
    
    % save data
    evalsList{ss} = evals;
    hvalsList{ss} = hvals;
    
    % voronoi
    fprintf('* Computing Voronoi data for sim ...\n');
    [a, ~, voroAreas, voroCalA] = getSurfaceVoronoi(xpos,ypos,l0,nv,L(1));
    fprintf('...Voronoi done!\n');
    
    aList{ss} = a;
    voroAreasList{ss} = voroAreas;
    voroCalAList{ss} = voroCalA;
end

% remove problematic sims
simList(simSkip)        = [];
NCELLSList(simSkip)     = [];
NvList(simSkip)         = [];
LList(simSkip,:)        = [];
zcList(simSkip)         = [];
zvList(simSkip)         = [];
a0List(simSkip)         = [];
l0List(simSkip)         = [];
calAList(simSkip,:)     = [];
aList(simSkip,:)        = [];
voroAreasList(simSkip,:) = [];
voroCalAList(simSkip,:) = [];
phiJList(simSkip,:)     = [];
evalsList(simSkip)      = [];
hvalsList(simSkip)      = [];
pratioList(simSkip)     = [];

% save to matfile
save(saveStr,'simList','NCELLSList','NvList','LList','zcList','zvList','a0List',...
    'l0List','calAList','phiJList','evalsList','hvalsList','pratioList',...
    'aList','voroAreasList','voroCalAList');


end


function [a, calA, voroAreas, voroCalA] = getSurfaceVoronoi(xpos,ypos,l0,nv,L) 

%% Get Voronoi diagram from Delaunay triangulation

% number of particles
NCELLS = length(nv);

% all vertices + interpolated points
NINTERP = 15;
NVTOT = (NINTERP+1)*sum(nv);

xall = zeros(sum(nv),1);
yall = zeros(sum(nv),1);
gi = 1;
for cc = 1:NCELLS
    xtmp = xpos{cc};
    ytmp = ypos{cc};
    for vv = 1:nv(cc)
        xall(gi) = xtmp(vv);
        yall(gi) = ytmp(vv);
        gi = gi + 1;
    end
end

% get vertices in global coordinates
gx = zeros(9*NVTOT,1);
gy = zeros(9*NVTOT,1);
ci = zeros(NVTOT,1);
gi = 1;
for nn = 1:NCELLS
    ip1 = [2:nv(nn) 1];
    xtmp = xpos{nn};
    ytmp = ypos{nn};
    for vv = 1:nv(nn)
        % add main point
        gx(gi) = xtmp(vv);
        gy(gi) = ytmp(vv);
        ci(gi) = nn;
        
        if NINTERP == 0
            gi = gi + 1;
        else
            % get segment to next
            lx = xtmp(ip1(vv)) - xtmp(vv);
            ly = ytmp(ip1(vv)) - ytmp(vv);
            l = sqrt(lx^2 + ly^2);
            del = l/(NINTERP+1);
            ux = lx/l;
            uy = ly/l;
            gi = gi + 1;
            for ii = 1:NINTERP
                gx(gi) = gx(gi-1) + del*ux;
                gy(gi) = gy(gi-1) + del*uy;
                ci(gi) = nn;
                gi = gi + 1;
            end
        end
    end
end

% also add periodic boundaries
blk = 1;
for xx = -1:1
    for yy = -1:1
       if (xx == 0 && yy == 0)
           continue;
       end
       
       % indices
       i0 = NVTOT*blk + 1;
       i1 = NVTOT*(blk + 1);
       
       % add
       gx(i0:i1) = gx(1:NVTOT) + xx*L;
       gy(i0:i1) = gy(1:NVTOT) + yy*L;
       
       % update block
       blk = blk + 1;
    end
end

% make delaunay
DT = delaunayTriangulation(gx,gy);

% voronoi diagram
[V,e] = voronoiDiagram(DT);

voroAreas = zeros(NCELLS,1);
svoroEdgeInfo = cell(NCELLS,1);
for vv = 1:NVTOT
    vvi = e{vv};
    vinfo = V(vvi,:);
    
    % add to area
    atmp = polyarea(vinfo(:,1),vinfo(:,2));
    voroAreas(ci(vv)) = atmp + voroAreas(ci(vv));
    
    % label all vertices as bndry, remove by checking interior edges
    NVCE = length(vvi);
    onbound = true(NVCE,1);
    civ = ci(vv);
    
    % find exterior surface by checking neighbors    
    vp1 = [vvi(2:end) vvi(1)];
    for ee = 1:NVCE
        % get pairs of vertices on edge
        ve1 = vvi(ee);
        ve2 = vp1(ee);
        
        % loop over other voronoi cells in this particle
        vint = find(ci == civ);
        NVINT = length(vint);
        for vvv = 1:NVINT
            if vint(vvv) == vv
                continue;
            end
            vvj = e{vint(vvv)};
            if sum(ve1 == vvj) == 1 && sum(ve2 == vvj) == 1
                onbound(ee) = false;
                break;
            end
        end
    end
    
    % add to master list
    for ee = 1:NVCE
        if onbound(ee)
            ve1 = vvi(ee);
            ve2 = vp1(ee);
            
            % add to face list
            svoroEdgeInfo{ci(vv)} = [svoroEdgeInfo{ci(vv)}; ve1, ve2];
        end
    end
end

% clean face list
svoroFaceList = cell(NCELLS,1);
for cc = 1:NCELLS
    % face list
    etmp = svoroEdgeInfo{cc};
    NVE = size(etmp,1);
    
    % construct face list
    ftmp = zeros(NVE,1);
    ftmp(1) = etmp(1,1);
    ftmp(2) = etmp(1,2);
    curr = 1;
    for ee = 2:NVE
        nxtind = etmp(curr,2) == etmp(:,1);
        ftmp(ee) = etmp(nxtind,1);
        curr = find(nxtind);
    end
    
    % save 
    svoroFaceList{cc} = ftmp;
end

% get particle areas
a = zeros(NCELLS,2);
calA = zeros(NCELLS,1);
voroCalA = zeros(NCELLS,1);
for cc = 1:NCELLS
    % cell shape parameters
    xtmp = xpos{cc};
    ytmp = ypos{cc};
    atmp = polyarea(xpos{cc},ypos{cc});
    ip1 = [2:nv(cc) 1];
    lx = xtmp(ip1) - xtmp;
    ly = ytmp(ip1) - ytmp;
    l = sqrt(lx.^2 + ly.^2);
    ptmp = sum(l);
    
    a(cc,1) = atmp;
    a(cc,2) = atmp + 0.25*pi*(l0(cc)^2)*(0.5*nv(nn)-1);
    calAn = nv(cc)*tan(pi/nv(cc))/pi;
    calA(cc) = (ptmp^2/(4.0*pi*atmp))/calAn;
    
    % get info for cell
    ftmp = svoroFaceList{cc};
    vatmp = voroAreas(cc);
    lx = V([ftmp(2:end); ftmp(1)],1) - V(ftmp,1);
    ly = V([ftmp(2:end); ftmp(1)],2) - V(ftmp,2);
    l = sqrt(lx.^2  + ly.^2);
    vptmp = sum(l);
    
    voroCalA(cc) = vptmp^2/(4.0*pi*vatmp);
end


end


function V2_norm = ModeProj_DPM(Nc, Ns, Dc, x, y, eigV)
%% FUNCTION to compute projection of modes onto different directions
% -- Nc: number of cells
% -- Ns: number of vertices / cell
% -- Dc: effective cell diamters, based on l0
% -- x: all x vertex positions
% -- y: all y vertex positions
% -- eigV: eigenVector matrix, sorted by all x, then all y

% jack's edits to Dongs code
N = sum(Ns);

ift = zeros(N, 1);
jft = zeros(N, 1);
idx_start = zeros(Nc, 1);
idx_end = zeros(Nc, 1);
Ns_all = zeros(N, 1);

for nc = 1:Nc
    idx1 = sum(Ns(1:nc-1)) + 1;
    idx2 = sum(Ns(1:nc));
    idx_start(nc) = idx1;
    idx_end(nc) = idx2;
    ift(idx1:idx2) = [2:Ns(nc) 1]' + idx1 - 1;
    jft(idx1:idx2) = [Ns(nc) 1:Ns(nc)-1]' + idx1 - 1;
    Ns_all(idx1:idx2) = Ns(nc);
end

U = zeros(2 * N, 3 * Nc);
for it = 1:Nc
    idx1 = idx_start(it);
    idx2 = idx_end(it);
    U(idx1:idx2, it) = 1;
    U(idx1+N:idx2+N, it + Nc) = 1;
    
    dtheta = 2 / Dc(nc);
    xs = x(idx1:idx2);
    ys = y(idx1:idx2);
    x_cen = mean(xs);
    y_cen = mean(ys);
    rx = xs - x_cen;
    ry = ys - y_cen;
    r = sqrt(rx.^2 + ry.^2);
    theta_s = atan2(ry, rx);
    U(idx1:idx2, it + 2 * Nc) = -dtheta * r .* sin(theta_s);
    U(idx1+N:idx2+N, it + 2 * Nc) = dtheta * r .* cos(theta_s);
end

for it = 1:3*Nc
    U(:, it) = U(:, it) / norm(U(:, it));
end

V2 = zeros(3 * Nc, 2 * N);

for i = 1:2*N
    for j = 1:3*Nc
        V2(j, i) = dot(eigV(:, i), U(:, j));
    end
end

V2_norm = zeros(3, 2*N);
for it = 1:2*N
    V2_norm(1, it) = sum(V2(1:2*Nc, it).^2); % translation
    V2_norm(2, it) = sum(V2(2*Nc+1:3*Nc, it).^2); % rotation 
    V2_norm(3, it) = 1 - V2_norm(1, it) - V2_norm(2, it); % deformation
end



end
 