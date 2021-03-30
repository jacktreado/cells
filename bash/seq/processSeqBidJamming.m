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
phiJList            = cell(NSIM,2);         % packing fractions (both from a and a0) at jamming

% VDOS list
evalsList           = cell(NSIM,1);         % list of eigenvalues at jamming
hvalsList           = cell(NSIM,1);         % list of projections of eigenvectors onto H

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

    % close the file
    fclose(fid);
    
    % save data
    evalsList{ss} = evals;
    hvalsList{ss} = hvals;
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
phiJList(simSkip,:)     = [];
evalsList(simSkip)      = [];
hvalsList(simSkip)      = [];

% save to matfile
save(saveStr,'simList','NCELLSList','NvList','LList','zcList','zvList','a0List',...
    'l0List','calAList','phiJList','evalsList','hvalsList');


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