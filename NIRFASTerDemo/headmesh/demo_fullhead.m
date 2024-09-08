clear;
% Here loaded is a segmented FEM head mesh
% no optodes are specified yet, so NIRFAST will throw a warning
mesh = load_mesh('Example_Mesh');
%% Let's set the optical properties for each tissue type
% These are for 850nm; from paper https://pubmed.ncbi.nlm.nih.gov/22330315/
% You need to use different numbers of different wavelengths
% mua in mm-1 and musp in mm-1
% csf
idx = mesh.region==1;
mesh.mua(idx) = 0.0040;
mesh.mus(idx) = 0.3;
% white
idx = mesh.region==2;
mesh.mua(idx) = 0.0208;
mesh.mus(idx) = 1.0107;
% gray
idx = mesh.region==3;
mesh.mua(idx) = 0.0192;
mesh.mus(idx) = 0.6726;
% bone
idx = mesh.region==4;
mesh.mua(idx) = 0.0139;
mesh.mus(idx) = 0.84;
% skin
idx = mesh.region==5;
mesh.mua(idx) = 0.0190;
mesh.mus(idx) = 0.64;

mesh.ri = 1.4*ones(size(mesh.ri));

% Define the optodes array
optodes = [
    [-75.497, -60.6597, 23.3938];
    [-76.527, -59.7415, 1.8209];
    [-79.543, -50.3532, 2.309];
    [-78.713, -49.7402, 23.3223];
    [-82.487, -40.9129, 1.9661];
    [-83.493, -30.3324, 1.1809];
    [-82.615, -20.173, 1.3609];
    [-81.516, -10.394, 1.599];
    [-79.482, -1.479, 0.1183];
    [-81.556, -39.6442, 23.9549];
    [-82.46, -30.7959, 22.5415];
    [-81.458, -20.068, 22.293];
    [-80.274, -9.984, 21.2728];
    [-77.486, 0.937, 20.5718];
    [-79.49, -50.6068, 12.0683];
    [-81.452, -40.0453, 11.1484];
    [-83.534, -29.9768, 10.9893];
    [-82.522, -20.54, 10.9675];
    [-81.297, -10.417, 10.5064];
    [-78.752, -0.313, 11.8946];
    [-76.728, -57.6585, 11.356];
    [-72.467, -70.0335, 2.1359];
    [-67.265, -80.3862, 1.9196];
    [-72.416, -70.0897, 12.3359];
    [-66.53, -80.3672, 11.6999];
    [-71.202, -70.3243, 23.2168];
    [-65.282, -80.2498, 22.9951];
    [-78.471, -0.145, -9.789];
    [-80.723, -10.806, -9.7021];
    [-82.44, -19.627, -9.697];
    [-83.537, -29.6207, -10.1978];
    [-81.784, -40.6013, -10.2447];
    [-78.927, -51.2362, -10.4246];
    [-75.544, -60.8856, -10.0527];
    [-72.445, -71.0607, -9.3775];
    [-65.961, -82.5101, -9.6401]
];

% Number of source and detectors 
nOptodes = size(optodes, 1);

%source structure
srcpos = optodes;
source = [];
source.coord = optodes; %optodes as sources
source.fixed = 0;
source.fwhm = zeros(nOptodes, 1);
source.num = (1:nOptodes)';


%detector structure
detpos = optodes;
detector = [];
detector.coord = optodes;
detector.fixed = 0;
detector.fwhm = zeros(nOptodes, 1);
detector.num = (1:nOptodes)';

%generating link matrix(creating all possible link excluding self link) 
[src_idx,det_idx] = ndgrid(1:nOptodes, 1:nOptodes);
valid_pairs = src_idx ~= det_idx;
link = [src_idx(valid_pairs), det_idx(valid_pairs), ones(sum(valid_pairs(:)),1)];

%mesh structure
mesh.source = source;
mesh.meas = detector;
mesh.link = link;


save_mesh('headmodel', mesh);

mesh = load_mesh('headmodel');

xgrid = -88:2:88;
ygrid = -118:2:84;
zgrid = -74:2:100;
mesh = gen_intmat(mesh, xgrid, ygrid, zgrid);
%jacobian on the grid
[J, data0] = jacobiangrid_stnd_FD(mesh);


% Have a look
figure, plot3dmesh(mesh);
hold on
scatter3(srcpos(:,1),srcpos(:,2),srcpos(:,3), 50, 'ro', 'filled')
scatter3(detpos(:,1),detpos(:,2),detpos(:,3), 50, 'bs', 'filled')

%% Add some activation and reconstruct it

% An area of activation and assume gray matter nodes within a radius
% Defining the starting point and dimensions
start_point = [-68.9658, -32.6004, 17.4106];
length_x = 20; 
length_y = 10; 
%depth from the brain surface
depth_from_surface = 10; 

% Defining the ranges for x and y coordinates based on the start point
x_range = [start_point(1), start_point(1) + length_x];
y_range = [start_point(2), start_point(2) + length_y];
z_value = start_point(3) - depth_from_surface; % The depth from the brain surface

% Finding nodes within the rectangular
in_x_range = mesh.nodes(:,1) >= x_range(1) & mesh.nodes(:,1) <= x_range(2);
in_y_range = mesh.nodes(:,2) >= y_range(1) & mesh.nodes(:,2) <= y_range(2);
near_z_value = abs(mesh.nodes(:,3) - z_value) <= 1;

% Combining conditions to find nodes
idx = in_x_range & in_y_range & near_z_value;

% Calculate the center of the activation region
activation_coords = mesh.nodes(idx, :);
activation_center = mean(activation_coords, 1);

% Updating the mesh properties for particular nodes
mesh2 = mesh;
mesh2.mua(idx) = 1.1 * mesh.mua(idx);

% Simulate measured data
data = femdata_stnd_FD(mesh2, 0);
% Reconstruction using Tikhonov regularization
delta_mua_recon = tikhonov(J.complete, .01, log(data.amplitude./data0.amplitude)+randn(size(data.amplitude))*0.5);

% If you have NeuroDOT and it's added to path, we can visualize the result
% very nicely by interpolating it onto a cortical surface mesh
[~,infoB]=LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti',[],'nii');
load('MNI164k_big.mat')
% The struct below is related to the original head model and the grid
% Calculation of it is non-intuitive so let's take the numbers for now
infoA = [];
infoA.center = [-91,-87,-103];
infoA.nVx = 89;
infoA.nVy = 102;
infoA.nVz = 88;
infoA.mmppix = [-2, -2, -2];

% Convert ground truth to grid
delta_mua0 = mesh.vol.mesh2grid*(mesh2.mua-mesh.mua);
delta_mua0 = reshape(delta_mua0, length(ygrid), length(xgrid), length(zgrid));
delta_mua0 = permute(delta_mua0, [2,1,3]);
% Interpolate to atlas
truth_atlas = affine3d_img(delta_mua0,infoA,infoB,eye(4));
pS.Scale=0.8*max(abs(truth_atlas(:)));     % Scale wrt/max of data
pS.Th.P=0;            % Threshold to see strong activations
pS.Th.N=0; 
pS.view='post'; % Posterior view
pS.ctx='std'; % Standard pial cortical view
PlotInterpSurfMesh(truth_atlas, MNIl,MNIr, infoB, pS);
%hold on
%scatter3(srcpos(:,1),srcpos(:,2),srcpos(:,3), 50, 'ro', 'filled')
%scatter3(detpos(:,1),detpos(:,2),detpos(:,3), 50, 'bs', 'filled')



% do the same thing with the reconstructed result
delta_mua_recon = reshape(delta_mua_recon, length(ygrid), length(xgrid), length(zgrid));
delta_mua_recon = permute(delta_mua_recon, [2,1,3]);
% Interpolate to atlas
recon_atlas = affine3d_img(delta_mua_recon,infoA,infoB,eye(4));
pS.Scale=0.8*max(abs(recon_atlas(:)));     % Scale wrt/max of data
pS.Th.P=0;            % Threshold to see strong activations
pS.Th.N=0; 
pS.view='post'; % Posterior view
pS.ctx='std'; % Standard pial cortical view
PlotInterpSurfMesh(recon_atlas, MNIl,MNIr, infoB, pS);
%% Calculating sensitivity of each cahnnel 
sensitivity = sum(abs(J.complete(:, idx)), 2);

% Normalize sensitivity
sensitivity = sensitivity / max(sensitivity);

%  histogram of the ROI sensitivity
figure;
histogram(sensitivity, 30); % 30 bins, adjust the number if needed

threshold = 0.5;

%gathering the high sensitivity channel based on given threshold 

high_sensitivity_indices = find(sensitivity > threshold);

src_coord = [];
det_coord = [];
% Sort channels by sensitivity
%[sorted_sensitivity, sorted_indices] = sort(high_sensitivity_indices, 'descend');


% Display sorted channel numbers along with their source-detector pairs and positions
fprintf('Channel\tSource Index\tDetector Index\tSource Position\t\t\tDetector Position\t\t\tSensitivity\n');
for i = 1:length(high_sensitivity_indices)
    channel_num = high_sensitivity_indices(i);
    src_num = link(channel_num, 1);
    det_num = link(channel_num, 2);
    fprintf('%d\t%d\t\t%d\t\t[%.2f, %.2f, %.2f]\t\t[%.2f, %.2f, %.2f]\t\t%.4f\n', ...
            channel_num, src_num, det_num, ...
            srcpos(src_num, 1), srcpos(src_num, 2), srcpos(src_num, 3), ...
            detpos(det_num, 1), detpos(det_num, 2), detpos(det_num, 3), ...
            sensitivity(channel_num));
    src_coord = [src_coord; srcpos(src_num, :)];
    det_coord = [det_coord; detpos(det_num, :)];
end

% Plot the head mesh
figure, plot3dmesh(mesh);
hold on;

% Plot the source optodes for high sensitivity channels
scatter3(src_coord(:,1), src_coord(:,2), src_coord(:,3), 50, 'ro', 'filled');
scatter3(det_coord(:,1), det_coord(:,2), det_coord(:,3), 50, 'y', 'filled');
% Plot the head mesh
figure, plot3dmesh(mesh);
hold on;
% Plot the detector optodes for high sensitivity channels
scatter3(det_coord(:,1), det_coord(:,2), det_coord(:,3), 50, 'y', 'filled');

title('High Sensitivity Source-Detector Pairs');
hold off;









% Plot only the source-detector pairs with high sensitivity using plot3
for i = 1:length(high_sensitivity_indices)
    channel_num = high_sensitivity_indices(i);
    src_num = link(channel_num, 1);
    det_num = link(channel_num, 2);
    
    % Get the coordinates for the source and detector
    %src_coord = srcpos(src_num, :);
    %det_coord = detpos(det_num, :);
    
    %figure, plot3dmesh(mesh);
    %hold on

    % Plot a line between the source and detector
    %plot3([src_coord(1), det_coord(1)], [src_coord(2), det_coord(2)], [src_coord(3), det_coord(3)], 'g-', 'LineWidth', 2);
    
    % Plot the source and detector as well (if desired)
    %scatter3(src_coord(1), src_coord(2), src_coord(3), 50, 'ro', 'filled'); % Source position
    %scatter3(det_coord(1), det_coord(2), det_coord(3), 50, 'bs', 'filled'); % Detector position
end


figure, plot3dmesh(mesh);
hold on

scatter3(src_coord(1), src_coord(2), src_coord(3), 50, 'ro', 'filled'); % Source position
scatter3(det_coord(1), det_coord(2), det_coord(3), 50, 'bs', 'filled'); % Detector position


%title('Source-Detector Pairs with Sensitivity > 0.5');
hold off;