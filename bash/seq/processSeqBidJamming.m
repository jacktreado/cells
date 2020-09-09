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
hProjList           = cell(NSIM,1);         % list of projections of eigenvectors onto H
sProjList           = cell(NSIM,1);         % list of projections of eignevectors onto S
% dirProjList         = cell(NSIM,3);         % list of projections onto T, R, S directions
% vertPrList          = cell(NSIM,1);         % vertex-based participation ratio at jamming
% cellPrList          = cell(NSIM,1);         % cell-based participation ratio at jamming


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
    end
    phiJ = sum(a)/(L(1)*L(2));
    phi0J = sum(a0)/(L(1)*L(2));
    
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
    if isempty(evals)
        fprintf('Evals not read, skipping\n');
        simSkip(ss) = true;
        fclose(fid);
        continue;
    else
        evals = evals{1};
    end
    
    % TMP: stiffness eigenvalues
    hvals = textscan(fid,'%f',dof);
    if isempty(hvals)
        fprintf('hvals not read, skipping\n');
        simSkip(ss) = true;
        fclose(fid);
        continue;
    else
        hvals = hvals{1};
    end
    

    % eigenvector data
    frmt = repmat('%f ',1,dof);
    evecs = textscan(fid,frmt,dof);
    if isempty(evecs)
        fprintf('Evecs not read, skipping\n');
        simSkip(ss) = true;
        fclose(fid);
        continue;
    else
        evecs = cell2mat(evecs);
    end

%     TMP
%     hproj = zeros(dof,1);
%     sproj = zeros(dof,1);
%     for dd = 1:dof
%         if feof(fid) == 1
%             fprintf('EOF found before end of projections, skipping...\n');
%             simSkip(ss) = true;
%             break;
%         else
%             data = textscan(fid,'%f',1);
%             if isempty(data{1})
%                 fprintf('Projections not read, skipping...\n');
%                 simSkip(ss) = true;
%                 break;
%             else
%                 hproj(dd) = data{1};
%             end
%         end
% 
%         if feof(fid) == 1
%             fprintf('EOF found before end of projections, skipping...\n');
%             simSkip(ss) = true;
%             break;
%         else
%             data = textscan(fid,'%f',1);
%             if isempty(data{1})
%                 fprintf('Projections not read, skipping...\n');
%                 simSkip(ss) = true;
%                 break;
%             else
%                 sproj(dd) = data{1};
%             end
%         end
%     end
%     if simSkip(ss) == true
%         fclose(fid);
%         continue;
%     end
%     TMP

    % close the file
    fclose(fid);
    
    % save data
    evalsList{ss} = evals;
    hProjList{ss} = hvals;
%     hProjList{ss} = hproj;
%     sProjList{ss} = sproj;
    
    
    % NOTE IF YOU WANT TO SAVE PARTICIPATION RATIOS, PROJECTION ONTO T, R, etc
    % DO so here using VDOS data
end

% remove problematic sims
simList(simSkip) = [];
NCELLSList(simSkip) = [];
NvList(simSkip) = [];
LList(simSkip,:) = [];
zcList(simSkip) = [];
zvList(simSkip) = [];
a0List(simSkip) = [];
l0List(simSkip) = [];
calAList(simSkip,:) = [];
phiJList(simSkip,:) = [];
evalsList(simSkip) = [];
hProjList(simSkip) = [];
sProjList(simSkip) = [];

% save to matfile
save(saveStr,'simList','NCELLSList','NvList','LList','zcList','zvList','a0List',...
    'l0List','calAList','phiJList','evalsList','hProjList','sProjList');


end