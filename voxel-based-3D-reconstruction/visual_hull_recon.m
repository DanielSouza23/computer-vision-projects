clc; 
clear;

x_min = -2.5;
x_max =  2.5;
y_min = -3.0;
y_max =  3.0;
z_min =  0.0;
z_max =  2.5;

% grid resolution
Nx = 100;
Ny = 100;
Nz = 100;

% Spacing
x_coords = linspace(x_min, x_max, Nx);
y_coords = linspace(y_min, y_max, Ny);
z_coords = linspace(z_min, z_max, Nz);

% 3D grid
[X, Y, Z] = ndgrid(x_coords, y_coords, z_coords);

voxel_occupancy = zeros(Nx, Ny, Nz, 'uint8');

% Camera calibration 
rawP = [ 776.649963  -298.408539 -32.048386  993.1581875 132.852554  120.885834  -759.210876 1982.174000 0.744869  0.662592  -0.078377 4.629312012;
    431.503540  586.251892  -137.094040 1982.053375 23.799522   1.964373    -657.832764 1725.253500 -0.321776 0.869462  -0.374826 5.538025391;
    -153.607925 722.067139  -127.204468 2182.4950   141.564346  74.195686   -637.070984 1551.185125 -0.769772 0.354474  -0.530847 4.737782227;
    -823.909119 55.557896   -82.577644  2498.20825  -31.429972  42.725830   -777.534546 2083.363250 -0.484634 -0.807611 -0.335998 4.934550781;
    -715.434998 -351.073730 -147.460815 1978.534875 29.429260   -2.156084   -779.121704 2028.892750 0.030776  -0.941587 -0.335361 4.141203125;
    -417.221649 -700.318726 -27.361042  1599.565000 111.925537  -169.101776 -752.020142 1982.983750 0.542421  -0.837170 -0.070180 3.929336426;
    94.934860   -668.213623 -331.895508 769.8633125 -549.403137 -58.174614  -342.555359 1286.971000 0.196630  -0.136065 -0.970991 3.574729736;
    452.159027  -658.943909 -279.703522 883.495000  -262.442566 1.231108    -751.532349 1884.149625 0.776201  0.215114  -0.592653 4.235517090];

P = zeros(3,4,8);

for i = 1:8
    for j = 1:3
        P(j,1:4,i) = rawP(i, 4*(j-1)+1 : 4*(j-1)+4);
    end
end

silhouettes = cell(8, 1);

for c = 0:7
    silhName = sprintf('silh_cam%02d_00023_0000008550.pbm', c);
    I = imread(silhName);      % binary
    silhouettes{c+1} = I > 0;  % silhouettes 1-8 = true/false for cam 0-7
end

for ix = 1:Nx
    for iy = 1:Ny
        for iz = 1:Nz
            % 3D voxel center in homogeneous coordinates
            A = [ X(ix,iy,iz); Y(ix,iy,iz); Z(ix,iy,iz); 1];
            vox = true; % intially assumes the voxel is inside the hull

            % Check voxel against all cameras
            for cam = 1:8
                Pcam = P(:,:,cam);   % 3x4
                % Project: X = Px
                x = Pcam*A;
                u = x(1) / x(3); % homogenize
                v = x(2) / x(3);

                % pixel coordinates
                u_px = round(u);
                v_px = round(v);

                [H, W] = size(silhouettes{cam});

                % If out of bounds or outside silh: not in visual hull
                if u_px < 1 || u_px > W || v_px < 1 || v_px > H || ~silhouettes{cam}(v_px, u_px)
                    vox = false;
                    break;
                end
            end

            % If voxel is inside silhouette in every image, occupied
            if vox
                voxel_occupancy(ix,iy,iz) = 1;
            end
        end
    end
end


% Surface Voxels
surface = zeros(Nx, Ny, Nz, 'uint8');

% avoid out-of-bounds neighbors
for ix = 2:Nx-1
    for iy = 2:Ny-1
        for iz = 2:Nz-1
            if voxel_occupancy(ix,iy,iz) == 1  % only check occupied voxels
                % 6 neighbors, 1 = occupied and 0 = empty - needs at least 1 zero
                n1 = voxel_occupancy(ix+1, iy, iz);
                n2 = voxel_occupancy(ix-1, iy, iz);
                n3 = voxel_occupancy(ix, iy+1, iz);
                n4 = voxel_occupancy(ix, iy-1, iz);
                n5 = voxel_occupancy(ix, iy, iz+1);
                n6 = voxel_occupancy(ix, iy, iz-1);
                if (n1 == 0) || (n2 == 0) || (n3 == 0) || (n4 == 0) || (n5 == 0) || (n6 == 0)
                    surface(ix,iy,iz) = 1;
                end
            end
        end
    end
end

% Spacing - based on linspace
dx = x_coords(2) - x_coords(1);
dy = y_coords(2) - y_coords(1);
dz = z_coords(2) - z_coords(1);

points = [];
voxel_idx = [];

for ix = 2:Nx-1
    for iy = 2:Ny-1
        for iz = 2:Nz-1

            if surface(ix,iy,iz) == 1 % if a surface voxel
                xc = X(ix,iy,iz);
                yc = Y(ix,iy,iz);
                zc = Z(ix,iy,iz);

                % if neighbor is EMPTY, add face center
                if voxel_occupancy(ix+1,iy,iz) == 0
                    points = [points; xc + dx/2, yc, zc];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
                if voxel_occupancy(ix-1,iy,iz) == 0
                    points = [points; xc - dx/2, yc, zc];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
                if voxel_occupancy(ix,iy+1,iz) == 0
                    points = [points; xc, yc + dy/2, zc];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
                if voxel_occupancy(ix,iy-1,iz) == 0
                    points = [points; xc, yc - dy/2, zc];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
                if voxel_occupancy(ix,iy,iz+1) == 0
                    points = [points; xc, yc, zc + dz/2];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
                if voxel_occupancy(ix,iy,iz-1) == 0
                    points = [points; xc, yc, zc - dz/2];
                    voxel_idx = [voxel_idx; ix, iy, iz];
                end
            end
        end
    end
end

threeD = size(points, 1);

Xpts = points(:,1);
Ypts = points(:,2);
Zpts = points(:,3);

X_min = min(Xpts);
X_max = max(Xpts);
Y_min = min(Ypts);
Y_max = max(Ypts);
Z_min = min(Zpts);
Z_max = max(Zpts);

% Compute channels 

R = 255*(Xpts - X_min)./(X_max - X_min);
G = 255*(Ypts - Y_min)./(Y_max - Y_min);
B = 255*(Zpts - Z_min)./(Z_max - Z_min);

colors = uint8([R, G, B]);

% False-color RGB

ply = 'dancer.ply';

fid = fopen(ply, 'w');

% Header
fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'element vertex %d\n', threeD);
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'element face 0\n'); 
fprintf(fid, 'end_header\n');

% one vertex per line: x y z r g b
for i = 1:threeD
    fprintf(fid, '%f %f %f %d %d %d\n', points(i,1), points(i,2), points(i,3), colors(i,1), colors(i,2), colors(i,3));
end
fclose(fid);


for i=1:8
    tempP = P(:,:,i);
    C(1:4,i) = null(tempP);
    C(1:4,i)=C(1:4,i)/C(4,i);
end

% RGB for each camera
images = cell(8,1);
for c = 0:7
    img = sprintf('cam%02d_00023_0000008550.png', c);
    rgb = imread(img);  
    images{c+1} = rgb;
end

% Visibility
TR = zeros(Nx,Ny,Nz);
TG = zeros(Nx,Ny,Nz);
TB = zeros(Nx,Ny,Nz);

steps = 10; 

for ix = 2:Nx-1
    for iy = 2:Ny-1
        for iz = 2:Nz-1

            if surface(ix,iy,iz) == 1
                xc = X(ix,iy,iz);
                yc = Y(ix,iy,iz);
                zc = Z(ix,iy,iz);

                R_ = [];
                G_ = [];
                B_ = [];

                % check each camera
                for cam = 1:8
                    % camera center in world coords
                    CW = C(1:3, cam);

                    % direction from voxel to camera
                    dir = CW - [xc; yc; zc];

                    % visibility test
                    visible = true;
                    for n = 1:steps-1
                        t = n/steps;
                        pos = [xc; yc; zc] + t * dir;
                      
                        % Convert world coordinate to nearest voxel
                        tx = (pos(1) - x_min) / dx;   
                        ty = (pos(2) - y_min) / dy;  
                        tz = (pos(3) - z_min) / dz; 
                        
                        ix_s = round(tx) + 1; % +1 due to matlab indexing
                        iy_s = round(ty) + 1;
                        iz_s = round(tz) + 1;

                        if ix_s < 1 || ix_s > Nx || iy_s < 1 || iy_s > Ny || iz_s < 1 || iz_s > Nz
                            break;
                        end

                        % check for occlusion hits
                        if voxel_occupancy(ix_s,iy_s,iz_s) == 1 && ~(ix_s == ix && iy_s == iy && iz_s == iz)
                            visible = false;
                            break;
                        end
                    end

                    if ~visible % Skip camera if vision is occluded
                        continue;  
                    end

                    % If visible, project voxel center into this image
                    Pcam = P(:,:,cam);
                    h = Pcam * [xc; yc; zc; 1];
                    u = h(1)/h(3);
                    v = h(2)/h(3);
                    % Pixels
                    u_px = round(u);
                    v_px = round(v);

                    [H, W, ~] = size(images{cam});

                    if u_px < 1 || u_px > W || v_px < 1 || v_px > H
                        continue;   % outside image
                    end

                    pixel = images{cam}(v_px, u_px, :);   
                    pixel = double(pixel);              

                    R_(end+1) = pixel(1);
                    G_(end+1) = pixel(2);
                    B_(end+1) = pixel(3);
                end

                if ~isempty(R_) % If at least one camera captured color, take median values
                    TR(ix,iy,iz) = median(R_);
                    TG(ix,iy,iz) = median(G_);
                    TB(ix,iy,iz) = median(B_);
                end
            end
        end
    end
end

Tcolors = zeros(threeD, 3, 'uint8');

for k = 1:threeD
    ix = voxel_idx(k,1);
    iy = voxel_idx(k,2);
    iz = voxel_idx(k,3);

    Tcolors(k,1) = uint8(TR(ix,iy,iz));
    Tcolors(k,2) = uint8(TG(ix,iy,iz));
    Tcolors(k,3) = uint8(TB(ix,iy,iz));
end

% True-color RGB

ply = 'Color_dancer.ply';

fid = fopen(ply, 'w');

% Header
fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'element vertex %d\n', threeD);
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'element face 0\n'); 
fprintf(fid, 'end_header\n');

% one vertex per line: x y z r g b
for i = 1:threeD
    fprintf(fid, '%f %f %f %d %d %d\n', points(i,1), points(i,2), points(i,3), Tcolors(i,1), Tcolors(i,2), Tcolors(i,3));
end
fclose(fid);