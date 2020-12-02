% George Durrant

%% Femoral Implant Stress Distribution

clear
clc

A = imread('images/Femur','png');
A = rgb2gray(A);  % Convert the image to grayscale
A = double(A);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision

res = 0.78e-3;	% resolution in meters per pixel

% Find BMD weighted centroid 
sumy = 0;
sumz = 0;
total_sum = 0;

for i = 1:86
	for j = 1:86
            z = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y = (i-1)*res+res/2;	% y position of pixel center in m
            sumz = sumz+z*A(i,j);
            sumy = sumy+y*A(i,j);
            total_sum = total_sum + A(i,j);
	end
end


% Show image
figure(1)
imagesc(A),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

% Centroid coodinates
y_hat = sumy/total_sum;	% in m
z_hat = sumz/total_sum;	% in m
% fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
% fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

% Plot centroid axes on the image 
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])

%% Setup

% Define constants
My = 46;    % Nm
Mz = 28;    % Nm
Ey = 17e9;   % Pa, young bone
Eo = 13e9;   % Pa, old bone
pix_area = res*res; % area of each pixel, m^2

% Set up matrix of modulus values
BMD = zeros(86,86);
YoungE = A;
OldE = A;

for i = 1:86
	for j = 1:86
        BMD(i,j) = A(i,j)*2/255;
        % if BMD less than 1 it is empty space
        if BMD(i,j) < 1
            YoungE(i,j) = 0;
            OldE(i,j) = 0;
        else
            YoungE(i,j) = Ey;
            OldE(i,j) = Eo;
        end
    end
end

% load implant image
C = imread('images/Implant', 'png');
C = rgb2gray(C);
C = double(C);

% convert grayscale to modulus values
D = C;
for i = 1:15
    for j = 1:15
        if C(i,j) > 0 
            D(i,j) = 200e9;
        end
    end
end

% implant coordinates
ai = 50;
aj = 37;
bi = 37;
bj = 35;
ci = 41;
cj = 43;
di = 34;
dj = 51;

% create a matrix for each implant, overlapping E values on femur E values
ImplantAY = YoungE;
ImplantAO = OldE;
implant_i = 1;
for i = ai-7:ai+7
    implant_j = 1;
    for j = aj-7:aj+7
        if D(implant_i, implant_j) > 0
            ImplantAY(i,j) = D(implant_i, implant_j);
            ImplantAO(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantBY = YoungE;
ImplantBO = OldE;
implant_i = 1;
for i = bi-7:bi+7
    implant_j = 1;
    for j = bj-7:bj+7
        if D(implant_i, implant_j) > 0
            ImplantBY(i,j) = D(implant_i, implant_j);
            ImplantBO(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantCY = YoungE;
ImplantCO = OldE;
implant_i = 1;
for i = ci-7:ci+7
    implant_j = 1;
    for j = cj-7:cj+7
        if D(implant_i, implant_j) > 0
            ImplantCY(i,j) = D(implant_i, implant_j);
            ImplantCO(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantDY = YoungE;
ImplantDO = OldE;
implant_i = 1;
for i = di-7:di+7
    implant_j = 1;
    for j = dj-7:dj+7
        if D(implant_i, implant_j) > 0
            ImplantDY(i,j) = D(implant_i, implant_j);
            ImplantDO(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

%% Young Bone No Implant

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = YoungE(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = YoungE(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

% initialize matrices for stress (in Pa and MPa)
ox = zeros(86,86); % Pa
oxY = zeros(86,86); % MPa

% sum over each pixel to calculate stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = YoungE(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxY(i,j) = ox(i,j)/1000000;
    end
end

% find max stress and location
[M, I] = max(oxY(:));
[I_row, I_col] = ind2sub(size(oxY),I);
M
max_y = ((I_row-1)*res + res/2)*1000
max_z = ((I_col-1)*res + res/2)*1000

% centroid location
fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

% plot centroid (and stresses?)
figure(2)
imagesc(A),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Young Bone No Implant'); axis square;
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
plot(z_hat/res, y_hat/res, 'ro')

%% Young Bone Implant A

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAY(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAY(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxAY = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantAY(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxAY(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxAY(i,j);
        if ImplantAY(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxY(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Young Bone Implant B

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBY(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBY(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxBY = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantBY(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxBY(i,j) = ox(i,j)/1000000;
    end
end

max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxBY(i,j);
        if ImplantBY(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxY(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Young Bone Implant C

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCY(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCY(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxCY = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantCY(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxCY(i,j) = ox(i,j)/1000000;
    end
end

max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxCY(i,j);
        if ImplantCY(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxY(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Young Bone Implant D

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDY(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDY(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxDY = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantDY(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxDY(i,j) = ox(i,j)/1000000;
    end
end

max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxDY(i,j);
        if ImplantDY(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxY(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Old Bone No Implant

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = OldE(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = OldE(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

% initialize matrices for stress (in Pa and MPa)
ox = zeros(86,86); % Pa
oxO = zeros(86,86); % MPa

% sum over each pixel to calculate stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = OldE(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxO(i,j) = ox(i,j)/1000000;
    end
end

% find max stress and location
[M, I] = max(oxO(:));
[I_row, I_col] = ind2sub(size(oxO),I);
M
max_y = ((I_row-1)*res + res/2)*1000
max_z = ((I_col-1)*res + res/2)*1000

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

% plot centroid (and stresses?)
figure(2)
imagesc(A),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Old Bone No Implant'); axis square;
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
plot(z_hat/res, y_hat/res, 'ro')

%% Old Bone Implant A

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAO(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAO(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxAO = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantAO(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxAO(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxAO(i,j);
        if ImplantAO(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxO(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Old Bone Implant B

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBO(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBO(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxBO = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantBO(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxBO(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxBO(i,j);
        if ImplantBO(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxO(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Old Bone Implant C

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCO(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCO(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxCO = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantCO(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxCO(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxCO(i,j);
        if ImplantCO(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxO(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Old Bone Implant D

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDO(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDO(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxDO = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantDO(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxDO(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxDO(i,j);
        if ImplantDO(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxO(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone Setup

% make two functions using given values of BMD, E
f = @(b) 17/(2^b);
g = @(b) 1.107/(0.26^b);

% find the intersection point
b = fzero(@(b)f(b)-g(b),0);
a = f(b);

% find the modulus values for each pixel
HetE = zeros(86,86);
for i = 1:86
    for j = 1:86
        HetE(i,j) = a*(BMD(i,j)^b)*1000000000; % Pa
    end
end

% create a matrix for each implant, overlapping E values on femur E values
ImplantAH = HetE;
implant_i = 1;
for i = ai-7:ai+7
    implant_j = 1;
    for j = aj-7:aj+7
        if D(implant_i, implant_j) > 0
            ImplantAH(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantBH = HetE;
implant_i = 1;
for i = bi-7:bi+7
    implant_j = 1;
    for j = bj-7:bj+7
        if D(implant_i, implant_j) > 0
            ImplantBH(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantCH = HetE;
implant_i = 1;
for i = ci-7:ci+7
    implant_j = 1;
    for j = cj-7:cj+7
        if D(implant_i, implant_j) > 0
            ImplantCH(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

ImplantDH = HetE;
implant_i = 1;
for i = di-7:di+7
    implant_j = 1;
    for j = dj-7:dj+7
        if D(implant_i, implant_j) > 0
            ImplantDH(i,j) = D(implant_i, implant_j);
        end
        implant_j = implant_j+1;
    end
    implant_i = implant_i+1;
end

%% Heterogeneous Bone No Implant

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = HetE(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = HetE(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

% initialize matrices for stress (in Pa and MPa)
ox = zeros(86,86); % Pa
oxH = zeros(86,86); % MPa

% sum over each pixel to calculate stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = HetE(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxH(i,j) = ox(i,j)/1000000;
    end
end

% find max stress and location
[M, I] = max(oxH(:));
[I_row, I_col] = ind2sub(size(oxH),I);
M
max_y = ((I_row-1)*res + res/2)*1000
max_z = ((I_col-1)*res + res/2)*1000

% find min stress for NA
min_stress_upper = 100;
min_stress_lower = 100;
min_i_upper = 0;
min_j_upper = 0;
min_i_lower = 0;
min_j_lower = 0;

for i = 1:86
    for j = 1:86
        stress = oxH(i,j);
        if i < 44
            if abs(stress) < abs(min_stress_upper) && stress ~= 0
                min_stress_upper = stress;
                min_i_upper = i;
                min_j_upper = j;
            end
        else 
            if abs(stress) < abs(min_stress_lower) && stress ~= 0 
                min_stress_lower = stress;
                min_i_lower = i;
                min_j_lower = j;
            end
        end
    end
end

slope_H = (min_j_lower - min_j_upper)/(min_i_lower - min_i_upper)
inter_H = min_i_lower - 1/slope_H*min_j_lower

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone Implant A

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAH(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantAH(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxAH = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantAH(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxAH(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxAH(i,j);
        if ImplantAH(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxH(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone Implant B

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBH(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantBH(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxBH = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantBH(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxBH(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxBH(i,j);
        if ImplantBH(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxH(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone Implant C

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCH(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantCH(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxCH = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantCH(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxCH(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

for i = 1:86
    for j = 1:86
        stress = oxCH(i,j);
        if ImplantCH(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxH(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone Implant D

sum_y = 0;
sum_z = 0;
sum_den = 0;

% sum over each pixel for y_hat, z_hat
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDH(i,j);                 % Young's modulus of pixel
            sum_y = sum_y + Ei*pix_area*y_bar;
            sum_z = sum_z + Ei*pix_area*z_bar;
            sum_den = sum_den + Ei*pix_area;
	end
end

y_hat = sum_y/sum_den;	% in m
z_hat = sum_z/sum_den;	% in m

% now do polar moment of inertia using y_hat, z_hat
I_bar = res^4/12; % same for Iyy_bar, Izz_bar, Iyz_bar

Iyy_hat = 0;
Izz_hat = 0;
Iyz_hat = 0;

Iyy_star = 0;
Izz_star = 0;
Iyz_star = 0;

% sum over each pixel for Iyy_star, Izz_star, Iyz_star
for i = 1:86
	for j = 1:86
            z_bar = (j-1)*res+res/2;	% z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;	% y position of pixel center in m
            Ei = ImplantDH(i,j);                % Young's modulus of pixel
            
            Iyy_hat = I_bar + pix_area*(z_bar - z_hat)^2;
            Izz_hat = I_bar + pix_area*(y_bar - y_hat)^2;
            Iyz_hat = pix_area*(z_bar - z_hat)*(y_bar - y_hat);
            
            Iyy_star = Iyy_star + Ei*Iyy_hat;
            Izz_star = Izz_star + Ei*Izz_hat;
            Iyz_star = Iyz_star + Ei*Iyz_hat;
	end
end

ox = zeros(86,86); % stress matrix
oxDH = zeros(86,86);

% sum over each pixel for stress
for i = 1:86
    for j= 1:86
            z_bar = (j-1)*res+res/2;    % z position of pixel CENTER in m
            y_bar = (i-1)*res+res/2;    % z position of pixel center in m
            Ei = ImplantDH(i,j);                % Young's modulus of pixel
            s = z_bar - z_hat;
            t = y_bar - y_hat;
            
            ox(i,j) = Ei*((My*Izz_star + Mz*Iyz_star)*s - (My*Iyz_star + Mz*Iyy_star)*t) / (Iyy_star*Izz_star - Iyz_star^2);
            oxDH(i,j) = ox(i,j)/1000000;
    end
end

% find max stress in bone
max_stress = 0;
max_i = 0;
max_j = 0;

% find min stress for NA
min_stress_upper = 100;
min_stress_lower = 100;
min_i_upper = 0;
min_j_upper = 0;
min_i_lower = 0;
min_j_lower = 0;

for i = 1:86
    for j = 1:86
        stress = oxDH(i,j);
        if ImplantDH(i,j) < 200e9
            if stress > max_stress
                max_stress = stress;
                max_i = i;
                max_j = j;
            end
            if i < 44
                if abs(stress) < abs(min_stress_upper) && stress ~= 0
                    min_stress_upper = stress;
                    min_i_upper = i;
                    min_j_upper = j;
                end
            else 
                if abs(stress) < abs(min_stress_lower) && stress ~= 0 
                    min_stress_lower = stress;
                    min_i_lower = i;
                    min_j_lower = j;
                end
            end
        end
    end
end

max_stress
max_y = ((max_i-1)*res + res/2)*1000
max_z = ((max_j-1)*res + res/2)*1000

slope_D = (min_j_lower - min_j_upper)/(min_i_lower - min_i_upper)
inter_D = min_i_lower - 1/slope_D*min_j_lower

% find normalized stress to No Implant case at this location
norm_stress = max_stress/oxH(max_i, max_j)*100

fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

%% Heterogeneous Bone - Implant D contour

[Y,Z] = meshgrid(1:86);
contourf(Y, Z, oxDH, 30)
colorbar
caxis([-40 80])
set(gca, 'YDir', 'reverse')

hold on
NA_D = @(x) 1/slope_D*x + inter_D
fplot(NA_D, 'red')
xlim([15 70])
ylim([15 70])


%% Heterogeneous Bone - No Implant contour

contourf(Y, Z, oxH, 30)
hold on
colorbar
caxis([-40 80])
set(gca, 'YDir', 'reverse')

hold on
NA_H = @(x) 1/slope_H*x + inter_H
fplot(NA_H, 'red')
xlim([15 70])
ylim([15 70])

