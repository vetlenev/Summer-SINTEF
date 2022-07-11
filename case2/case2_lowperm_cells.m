%% 2D buoyancy-driven CO2 migration - low permeable cells
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui
ROOTDIR = strrep(ROOTDIR, '\', '/');

%% Define 2D grid
nx = 150; ny = 1; nz = 150; % 100 200 | 30 60
lx = 300; ly = 1; lz = 300; % 300 400 | 50 100
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);
x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);

%% Define rock and fluid objects
perm = repmat(100, [G.cells.num,1])*milli*darcy;
poro = 0.5;

%% Directories
n_lowperm_layers = 40;

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case2/plots_residual');
data_dir = strcat(ROOTDIR, '../summer_sintef/case2/data_residual');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_lowperm_layers);
dir_exists = mkdir(plot_base_dir);
dir_exists = mkdir(data_dir);
dir_exists = mkdir(plot_dir);

%% Define low-perm cells
seed = rng(1010);
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);

rand_x_start = linspace(perc_val(x, 0.0), perc_val(x, 0.9), n_lowperm_layers);
rand_x_start = rand_x_start(randperm(length(rand_x_start))); % random permute

rand_x_stop = rand_x_start + randi([round(perc_val(x, 0.15)), ...
                                    round(perc_val(x, 0.45))], ...
                                    [n_lowperm_layers, 1]).';
rand_x_stop = min(rand_x_stop, round(perc_val(x, 1))); % prevent layer going out of bounds

anticline_idx = 1:3:n_lowperm_layers;
x_anticline_start = rand_x_start(anticline_idx); %2:2:n_lowperm_layers
x_anticline_stop = rand_x_stop(anticline_idx); % odd layers curved

line_idx = 1:n_lowperm_layers;
line_idx(anticline_idx) = []; % indices for line layers
x_line_start = rand_x_start(line_idx); % 1:2:n_lowperm_layers
x_line_stop = rand_x_stop(line_idx); % even layers horizontal

z_start = linspace(perc_val(z, 0.1), perc_val(z, 0.9), n_lowperm_layers);
%z_line_start = z_start(1:2:n_lowperm_layers);
%z_anticline_start = z_start(2:2:n_lowperm_layers);

z_stop = z_start + perc_val(z, 0.02);
%z_line_stop = z_stop(1:2:n_lowperm_layers);
%z_anticline_stop = z_stop(2:2:n_lowperm_layers);


%% Generate layers
lowperm_cells = {};
all_lowperm_cells = [];

for i=line_idx
    lowperm_cells{i} = G.cells.indexMap(x > rand_x_start(i) & ...
                                    x < rand_x_stop(i) & ...
                                    z > z_start(i) & ...
                                    z < z_stop(i));
    perm(lowperm_cells{i}) = 10*milli*darcy;
    all_lowperm_cells = cat(1, all_lowperm_cells, lowperm_cells{i});
end

% Make anticline forms

for i=anticline_idx
    num_x_cells = numel(G.cells.indexMap(x > rand_x_start(i) ... 
                        & x < rand_x_stop(i) ...
                        & z == 1));
    theta = linspace(pi/4, 3*pi/4, num_x_cells);
    r = (rand_x_stop(i) - rand_x_start(i))/2;
    z_anticline_start = -r*sin(theta) + z_start(i) + r*cos(pi/4); % -r*sin to get anticline (since z positive downwards)
    z_anticline_stop = -r*sin(theta) + z_stop(i) + r*cos(pi/4); % + r*cos(pi/4) to shift cap to bottom of spherical cap
    
    anticline_cells = {};
    for j=1:num_x_cells
        anticline_cells{j} = G.cells.indexMap(x > (rand_x_stop(i) - rand_x_start(i))*(j-1)/num_x_cells + rand_x_start(i) & ...
                                    x < (rand_x_stop(i) - rand_x_start(i))*j/num_x_cells + rand_x_start(i) & ...
                                    z > z_anticline_start(j) & ...
                                    z < z_anticline_stop(j));
        perm(anticline_cells{j}) = 10*milli*darcy;
        all_lowperm_cells = cat(1, all_lowperm_cells, anticline_cells{j});
    end
                               
end

%% Compute rock+fluid objects
rock = makeRock(G, perm, poro);
T = computeTrans(G, rock, 'Verbose', true);

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3); % densities: [water, CO2]
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, G.cells.num)';
pc_vals = Case2Functions.Pc(dummy_Sw, 0, 100*barsa, 500*barsa, 2);

region_table = {[dummy_Sw, zeros(numel(dummy_Sw), 1)], [dummy_Sw, pc_vals]}; % container for pc values in each region
region_idx = {setdiff(G.cells.indexMap, all_lowperm_cells), all_lowperm_cells}; % region to interpolate (rest, lowperm)
%fluid_pc.pcOW = @(S, varargin) interpReg(region_table, S, reg_idx); % in both regions, interpolate over saturations S
% or with same function in entire domain
%fluid_pc.pcOW = Case2Functions.Pc(dummy_Sw, 0, 100*barsa, 500*barsa, 2);
% or using Leverett-J pc to account for differences in permeability

%% Horizontal well
rate = 30*meter^3/day; % 25 | 3
% Put well one cell above bottom to avoid it interacting with bottom BC
well_h = 3; % cell perforations in vertical direction
perforation_idx = G.cells.indexMap(z < max(z) & z > max(z)-well_h*lz/nz & x < lx/5);
W = addWell([], G, rock, perforation_idx, ...
            'Type', 'rate', 'Val', rate, ...
            'InnerProduct', 'ip_tpf', ...
            'Radius', 0.1, 'Dir', 'x', ...
            'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
            'Name', 'P1');

%W.lims.bhp = 1*barsa;

%% Plot grid
f1 = figure(1);
plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, perm, 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn); colorbar;
title('Permeability field (milli*darcy)');
axis equal tight
view([0, 0])
drawnow
hold off

saveas(f1, strcat(plot_dir, '/perm_field'), 'png');

%% Set up solver
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
disp(model)

%% Boundary conditions and schedule
bc = []; % no-flux as default
pz = fluid.rhoWS * norm(gravity) * uniquetol(z); % hydrostatic pressure for water (uniquetol necessary to avoid some duplicates by round-off)
bc = pside(bc, G, 'EAST', pz, 'sat', [1 0]);

% tot_time = 1600*day(); % 1600 days
% n_steps = 75; % 75
% dt = repmat(tot_time/n_steps, [n_steps, 1]);
dt = rampupTimesteps(4000*day(), 50*day(), 8);
n_steps = numel(dt);

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

%% Initial condition
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(G, 100*barsa, [1,0]);
t = 0;
p_mean = zeros(max(kk), n_steps);

lowperm_faces = zeros([size(all_lowperm_cells), 6]); % 6 faces in 3d cartesian coords
for i=1:size(all_lowperm_cells)
    icell_start = G.cells.facePos(all_lowperm_cells(i));
    icell_stop = G.cells.facePos(all_lowperm_cells(i)+1)-1;
    lowperm_faces(i,:) = G.cells.faces(icell_start:icell_stop, 1);
end

f2 = figure(2); % to hold saturations

plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar; caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f2, strcat(plot_dir, '/sat_0'), 'png');

% f3 = figure(3); % to plot pressure
% 
% plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
% plotCellData(G, state.pressure, 'EdgeColor', 'none');
% plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
% colorbar;
% title({'Pressure distribution (Pascal)' ['Time: ', formatTimeRange(t)]});
% axis equal tight
% view([0, 0])
% drawnow
% 
% saveas(f3, strcat(plot_dir, '/pres_0'), 'png');

%% Copmute solutions
[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

dt_plot = 5; % plot every 5th iteration
vol_ratios = ones(numel(states)+1, 1);
top_settled_cells = G.cells.indexMap(z < min(z(all_lowperm_cells)));

for i=1:numel(states)
    t = t+dt(i);
    S = states{i}.s(:,1);
    assert(max(S) < 1+eps && min(S) > -eps);
    
    p = reshape(states{i}.pressure, [max(ii), max(kk)]);
    p_mean(:,i) = mean(p, 1).';
    
    figure(2);
    plotCellData(G, S, 'EdgeColor', 'none');
    plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
    colormap(flipud(winter)); colorbar; caxis([0 1]);
    title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow, pause(0.2)
    
%     figure(3);
%     plotCellData(G, states{i}.pressure, 'EdgeColor', 'none');
%     plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
%     colorbar;
%     title({'Pressure distribution (Pascal)' ['Time: ', formatTimeRange(t)]});
%     axis equal tight
%     view([0, 0])
%     drawnow, pause(0.2)
    
    if ~mod(i, dt_plot)
       filename_f2 = sprintf(strcat(plot_dir, '/sat_%d'), i);
       %filename_f3 = sprintf(strcat(plot_dir, '/pres_%d'), i);
       saveas(f2, filename_f2, 'png');
       %saveas(f3, filename_f3, 'png');
    end
   
    %[tot_vol, lowperm_vol] = Case2Functions.Co2MassRatio(dims, gridsize, all_lowperm_cells, 1-S);
    [tot_vol, top_vol] = Case2Functions.Co2MassRatio(dims, gridsize, top_settled_cells, 1-S, rock, fluid);
    vol_ratio = 1 - top_vol/tot_vol; 
    vol_ratios(i+1) = vol_ratio;
end

vol_ratio_filename = sprintf(strcat(data_dir, '/vol_ratio_layers_%d.mat'), n_lowperm_layers);
save(vol_ratio_filename, 'vol_ratios'); % save to compare for different nr low-perm layers

%% Compare volume ratios for different nr low-perm layers
data_files = dir(strcat(data_dir, '/vol_ratio_layers_*.mat'));
run = {};
f4 = figure(4);

for i=1:numel(data_files)
   run{i} = load(strcat(data_files(i).folder, '\', data_files(i).name), 'vol_ratios');
   lab = regexp(data_files(i).name, 'layers_\d+', 'match');
   lab = strrep(lab{1}, '_', ': ');
   plot(run{i}.vol_ratios, 'o', 'DisplayName', lab);
   hold on;
end

legend('Location', 'northeast');
xlabel('Time step');
title('Ratio of total CO2 volume residing in low-permeable layers');
drawnow
saveas(f4, strcat(plot_base_dir, '/vol_ratio_', date), 'png');
hold off


%% TESTING
theta = linspace(pi/4, 3*pi/4, 100);
r = 3;
r0 = 5;
xx = linspace(0, 2, 100);
yy = r*sin(theta) + 1 - r;
legend('Location', 'east')
strcat(plot_base_dir, '_test_', date)
