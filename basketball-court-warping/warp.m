clear; 
clc; 

x  = [244 45; 419 69; 286 316; 0 193]; % Court's corners, tracked in Paint

xp = [1 1; 940 1; 940 500; 1 500];

H = part_a(x, xp);

IW = part_b(H);

HL = part_c(x);

IW2 = part_b(HL);

imwrite(IW, 'top_court_point.png');
imwrite(IW2, 'top_court_line.png');

figure;

subplot(1,2,1);
imshow(IW);
title('Court (Points)');

subplot(1,2,2);
imshow(IW2);
title('Court (Lines)');

function H = part_a(x, xp)
    N = 4; % # of points

    % Normalization
    [T, x_norm] = normalize(x, 488, 366); % original image
    [Tp, xp_norm] = normalize(xp, 940, 500); % new image

    % Notation: P: prime (new image)
    %           norm: x~, normalized point

    % Build A, (2Nx9)
    A = zeros(2*N, 9);
    for i = 1:N
        xi = x_norm(1,i) / x_norm(3,i); yi = x_norm(2,i) / x_norm(3,i); % pixel coordinates from the source image with w = 1
        ui = xp_norm(1,i) / xp_norm(3,i); vi = xp_norm(2,i) / xp_norm(3,i); % pixel coordinates in the new image with w = 1
        % (x,y,1) so divide by w (3rd element)
        
        % Slide 88 DLT A matrix, xi' = ui, yi' = vi, wi' = 1
        A(2*i-1, :) = [0, 0, 0, -xi, -yi,-1, vi*xi, vi*yi, vi]; % Each point gets 2 rows (in order) 
        A(2*i, :) = [xi, yi, 1, 0, 0, 0, -ui*xi, -ui*yi, -ui]; % (point 1 gets row 1 and 2, 2 gets 3 and 4), assemblying 8 rows for A
    end

    % Solve Ah = 0 via SVD (h = last column of V), slide 91
    [~,~,V] = svd(A, 0);
    h = V(:,end);   % 9x1 vector
    Hn = zeros(3,3); % Make it a 3x3 normalized H
    % 3x3 H matrix
    Hn(1,1) = h(1);  Hn(1,2) = h(2);  Hn(1,3) = h(3);
    Hn(2,1) = h(4);  Hn(2,2) = h(5);  Hn(2,3) = h(6);
    Hn(3,1) = h(7);  Hn(3,2) = h(8);  Hn(3,3) = h(9);

    % Denormalize h
    H = (Tp\Hn)*T;
    H = H ./ H(3,3); % Divide everything by H33 to normalize it through H33 = 1
end

function [T, x_norm] = normalize(x, w, h) % (slides 93–95)
    
    % Slide 93
    M = [w+h, 0, w/2; 0, w+h, h/2; 0, 0, 1];
    T = inv(M);

    xH = [x, ones(size(x,1),1)]';   % right now x is a 2x3 and needs to be a 3x3 to multiply by T
    x_norm = T*xH;                 % Transformation, slide 95
end

function IW = part_b(H)
    % height, width, color channels (3 = RGB) of original image
    h = 366;
    w = 488;
    C = 3;
    I = imread('basketball-court.ppm');
 
    % blank canvas for new image
    wp = 940;
    hp = 500;
    IW = zeros(hp, wp, C, 'uint8'); % uint8 = 0,255

    % Check every pixel
    for v = 1:hp           
        for u = 1:wp       
            g = [u; v; 1];
            f = H\g; %Inverse warping, slide 79,80
            x_ = f(1)/f(3); y_ = f(2)/f(3); % Normalize thru w = 1 
            if x_ >= 1 && x_ <= w-1 && y_ >= 1 && y_ <= h-1 % ok to interpolate, has neighbors on all 4 "corners" (not on border)
                % Find integer of partial pixel and its neighbor and distances, slide 81
                x1=floor(x_); x2=x1+1;
                y1=floor(y_); y2=y1+1; 
                a = x_ - x1;
                b = y_ - y1;
                
                % Bilinear Interpolation formula
                for k = 1:C % RGB values
                    f_i_j = double(I(y1, x1, k)); % matlab uses rows,columns so y,x
                    f_i1_j = double(I(y1, x2, k)); 
                    f_i_j1 = double(I(y2, x1, k)); 
                    f_i1_j1 = double(I(y2, x2, k)); 
                    
                    f_x_y = (1-a)*(1-b)*f_i_j + a*(1-b)*f_i1_j + a*b*f_i1_j1 + (1-a)*b*f_i_j1;

                    IW(v, u, k) = uint8(f_x_y); % Assign a color value to position (v,u) - new image
                end
            end
        end
    end
end

function HL = part_c(x)
    xH = [x, ones(size(x,1),1)]'; % making it a 3x3
    % cross product of two dots = line [a b c] -> ax + by + c = 0
    L = zeros(4,3);
    L(1,:) = cross(xH(:,1), xH(:,2))'; %top -> clockwise rotation
    L(2,:) = cross(xH(:,2), xH(:,3))'; 
    L(3,:) = cross(xH(:,3), xH(:,4))'; 
    L(4,:) = cross(xH(:,4), xH(:,1))'; 

    % Horizontal line y = 1 therefore x = 0 and c needs to be equal to -1 
    % (y + 0 -1 = 0)(top), vertical line x = 1 (left), etc
    Lp = [0 1 -1; 1 0 -940; 0 1 -500; 1 0 -1];
    
    % Linearize line parameters through w = 1
    for i = 1:4
        L(i,:)= L(i,:)/L(i,3);
        Lp(i,:) = Lp(i,:)/Lp(i,3);
    end

    % A matrix
    A = zeros(8, 9);
    for i = 1:4
        xL = L(i,1);  yL = L(i,2); % x and y components for original img  
        uL  = Lp(i,1); vL  = Lp(i,2); % x and y components for new img
        
        A(2*i-1, :) = [ -uL,  0,  uL*xL,  -vL,  0,  vL*xL,  -1,  0,  xL ];
        A(2*i, :) = [  0, -uL,  uL*yL,   0, -vL,  vL*yL,   0, -1,  yL ];
    end

    % Solve Ah = 0 via SVD (h = last column of V), slide 91
    [~,~,V] = svd(A, 0);
    h = V(:,end);   % 9x1 vector
    HL = zeros(3,3);
    % 3x3 H matrix
    HL(1,1) = h(1);  HL(1,2) = h(2);  HL(1,3) = h(3);
    HL(2,1) = h(4);  HL(2,2) = h(5);  HL(2,3) = h(6);
    HL(3,1) = h(7);  HL(3,2) = h(8);  HL(3,3) = h(9);

    HL = HL ./ HL(3,3); % Divide everything by H33 to normalize it through H33 = 1

    %https://www.researchgate.net/publication/220845575_Combining_Line_and_Point_Correspondences_for_Homography_Estimation
end
