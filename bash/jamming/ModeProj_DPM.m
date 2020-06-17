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