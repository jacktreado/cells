function singleParticleRelaxation(NV,Kl,Kb,calA0param,savestr)
%% FUNCTION to relax a single particle using FIRE + gradient descent, report final configuration and eigenvalue data

% initialize random number generator
rng('shuffle');

% relaxation parameters
Ftol        = 5e-15;
dt0         = 0.001;
itmax       = 2e8;
plotskip    = 20000;
phi0        = 0.1;
calA0       = calA0param*NV*tan(pi/NV)/pi;

% indexing variables
dof         = 2*NV;
xinds       = 1:2:dof-1;
yinds       = 2:2:dof;
ip1         = [2:NV 1];
im1         = [NV 1:NV-1];

%% Initialize particle

% initial coordinates
dth = 2.0*pi/NV;
th = 0:dth:(2.0*pi - dth);
x0 = cos(th)';
y0 = sin(th)';

% shape parameters
a0 = polyarea(x0,y0);
l0 = sqrt(4.0*pi*calA0*a0)/NV;

% rescale lengths
r0 = sqrt(a0);
a0 = a0/(r0*r0);
l0 = l0/r0;
x0 = x0./r0;
y0 = y0./r0;

% box length
L = sqrt(a0/phi0);

% % move to box center
x0 = x0 + 0.5*L;
y0 = y0 + 0.5*L;

% perturbation scale
p0 = 5e-1;

% perturb vertices
drx = randn(NV,1);
dry = randn(NV,1);
drnorm = sqrt(sum(drx.^2 + dry.^2,2));

drx = drx./drnorm;
dry = dry./drnorm;

x0 = x0 + p0*l0.*drx;
y0 = y0 + p0*l0.*dry;

x = x0;
y = y0;


%% RUN FIRE FROM SCRIPT, SAVE FIRE VARIABLES

% force parameters
rho0            = sqrt(a0);                 % units: length
fa              = 1/rho0;                   % units: inv length, because of grad a
fl              = Kl*(rho0/l0);             % units: dim. less
fb              = Kb*(rho0/(l0*l0));        % units: inv length, because of length in force expression

% initialize velocites to zero
vx = zeros(NV,1);
vy = zeros(NV,1);

% initialize forces
fx = zeros(NV,1);
fy = zeros(NV,1);

% FIRE VARIABLES (hard code in)
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dtmax       = 10*dt0;
dtmin       = 0.01*dt0;
dt          = dt0;
NNEGMAX     = 500;
NDELAY      = 20;
npPos       = 0;
npNeg       = 0;
alpha       = alpha0;

% force check
fcheck = 10*Ftol;

% USE FIRE to relax forces
it = 0;
t = 0.0;

% print to console
fprintf('** Relaxing shape for calA0 = %0.5g...\n',calA0);

% loop over FIRE protocol while fcheck and kcheck are above minimum
while(fcheck > Ftol && it < itmax)
    % update iterate
    it = it + 1;
    
    % Step 1. calculate P, fnorm, vnorm
    P = sum(fx.*vx) + sum(fy.*vy);
    
    % Step 2. adjust simulation based on net motion of system
    if P > 0
        % increase positive counter
        npPos = npPos + 1;
        
        % reset negative counter
        npNeg = 0;
        
        % alter simulation if enough positive steps have been taken
        if (npPos > NDELAY)
            % change time step
            if (dt*finc < dtmax)
                dt = dt*finc;
            end
        end
        
        % decrease alpha
        alpha = alpha*falpha;
    else
        % reset positive counter
        npPos = 0;
        
        % increase negative counter
        npNeg = npNeg + 1;
        
        % check for stuck simulation
        if (npNeg > NNEGMAX)
            fprintf('Simulation negative for too long, ending program here.\n')
            error('FIRE did not converge');
        end
        
        % decrease time step if past initial delay
        if (it > NDELAY)
            % decrease time step
            if (dt*fdec > dtmin)
                dt = dt*fdec;
            end
            
            % reset alpha
            alpha = alpha0;
        end
        
        % take a half step backwards
        x = x - 0.5*dt*vx;
        y = y - 0.5*dt*vy;
        
        % reset velocities to 0
        vx = zeros(NV,1);
        vy = zeros(NV,1);
    end
    
    % Step 1 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;
    
    vnorm = sqrt(sum(vx.*vx) + sum(vy.*vy));
    fnorm = sqrt(sum(fx.*fx) + sum(fy.*fy));

    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*(fx./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*(fy./fnorm)*vnorm;
    end
    
    % plot FIRE information 
    if mod(it,plotskip) == 0
        fprintf('\nOn FIRE step %d\n',it);
        fprintf('\t ** F = %0.5g\n',fcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** P = %0.5g\n',P);
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** vnorm = %0.5g\n',vnorm);
        fprintf('\t ** fnorm = %0.5g\n',fnorm);
    end
    
    % do first verlet update for vertices (assume unit mass)
    x = x + dt*vx;
    y = y + dt*vy;
    
    % update area
    a = polyarea(x, y);
    
    % * * * * * * * * * * * * * * * * * *
    % calculate forces based on positions
    % * * * * * * * * * * * * * * * * * *
    
    % reset forces
    fx = zeros(NV,1);
    fy = zeros(NV,1);
    
    % -- perimeter force
    
    % segment vectors
    lvx = x(ip1) - x;
    lvy = y(ip1) - y;
    
    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    
    % segment unit vectos
    ulvx = lvx./l;
    ulvy = lvy./l;
    
    % segment strain
    dli = (l./l0) - 1.0;
    dlim1 = (l(im1)./l0) - 1.0;
    
    % perimeter force at this iteration
    flx = fl*(dli.*ulvx - dlim1.*ulvx(im1));
    fly = fl*(dli.*ulvy - dlim1.*ulvy(im1));
    
    % add to total force
    fx = fx + flx;
    fy = fy + fly;
    
    
    % -- area force
    areaStrain = (a/a0) - 1.0;
    
    % area force at this iteration
    fax = fa*0.5*areaStrain.*(y(im1) - y(ip1));
    fay = fa*0.5*areaStrain.*(x(ip1) - x(im1));
    
    % add to total force
    fx = fx + fax;
    fy = fy + fay;
    
    % -- bending force
    
    % s vectors
    six = lvx - lvx(im1);
    siy = lvy - lvy(im1);
    
    % bending force at this iteration
    fbx = fb*(2.0*six - six(im1) - six(ip1));
    fby = fb*(2.0*siy - siy(im1) - siy(ip1));
    
    % add to force
    fx = fx + fbx;
    fy = fy + fby;
    
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;
    
    % update Fcheck
    fcheck = sqrt(sum(fx.^2 + fy.^2)/dof);
    
    % update total time spent
    t = t + dt;
end

%% Compute final eigenvalues

% project onto hessian
fprintf('** Shape is relaxed, computing eigenvalues...\n');
[Ha, Sa, Hl, Sl, Hb, Sb, Hbb, Sbb] = beltModelVDOS(x,y,a0,l0,1.0,Kl,Kb,0.0);

% compute eigenvalues for DM decomposition
H           = Ha + Hl + Hb + Hbb;
S           = -1.*(Sa + Sl + Sb + Sbb);
M           = H - S;
[V, ~]      = eig(M);
m           = zeros(dof,1);
h           = zeros(dof,1);

% constant translation modes
eX = zeros(dof,1);
eY = zeros(dof,1);

eX(xinds) = sqrt(1/NV)*ones(NV,1);
eY(yinds) = sqrt(1/NV)*ones(NV,1);

% tangential rotation vector
eR = zeros(dof,1);
eR(xinds) = -(y-mean(y));
eR(yinds) = x - mean(x);
eR = eR./norm(eR);

% use GS to orthoganalize modes to 3 known zero modes
vM = zeros(dof);
vM(:,1) = eX;
vM(:,2) = eY;
vM(:,3) = eR;

m(1) = eX'*M*eX;    h(1) = eX'*H*eX;
m(2) = eY'*M*eY;    h(2) = eY'*H*eX;
m(3) = eR'*M*eR;    h(3) = eR'*H*eR;
for mm = 4:dof
    % project modes from full D matrix
    vM(:,mm) = V(:,mm);
    for kk = 1:mm-1
        vM(:,mm) = vM(:,mm) - sum(V(:,mm).*vM(:,kk))*vM(:,kk);
    end

    % normalize mode
    vM(:,mm) = vM(:,mm)./norm(vM(:,mm));
    
    % compute eigenvalues
    m(mm) = vM(:,mm)'*M*vM(:,mm);
    h(mm) = vM(:,mm)'*H*vM(:,mm);
end


%% Save

fprintf('** Saving data to %s, ending.\n',savestr);
save(savestr,'NV','x','y','vM','m','h','H','S','a0','l0','Kl','Kb','it','t');
    

end
