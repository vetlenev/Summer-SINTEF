%% 2D buoyancy-driven CO2 migration - low permeable cells
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui
ROOTDIR = strrep(ROOTDIR, '\', '/');
seed = rng();
seed.Seed = 2135; % Must set here, otherwise not updated

%% Define 2D grid
nx = 100; ny = 1; nz = 50; % 200 100
lx = 1000; ly = 1; lz = 350; % 1000 350
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
n_lowperm_layers = 15;

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case3/plots_steady_state');
data_dir = strcat(ROOTDIR, '../summer_sintef/case3/data_steady_state');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_lowperm_layers);
dir_exists = mkdir(plot_base_dir);
dir_exists = mkdir(data_dir);
dir_exists = mkdir(plot_dir);

%% Define low-perm cells
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);

rand_x_start = linspace(perc_val(x, 0.0), perc_val(x, 0.9), n_lowperm_layers);
rand_x_start = rand_x_start(randperm(length(rand_x_start))); % random permute

rand_x_stop = rand_x_start + randi([round(perc_val(x, 0.1)), ...
                                    round(perc_val(x, 0.4))], ...
                                    [n_lowperm_layers, 1]).';
rand_x_stop = min(rand_x_stop, round(perc_val(x, 1))); % prevent layer going out of bounds

anticline_idx = 1:3:n_lowperm_layers; % every third layer is anticline
x_anticline_start = rand_x_start(anticline_idx); %2:2:n_lowperm_layers
x_anticline_stop = rand_x_stop(anticline_idx); % odd layers curved

line_idx = setdiff(1:n_lowperm_layers, anticline_idx); % indices for line layers
x_line_start = rand_x_start(line_idx); % 1:2:n_lowperm_layers
x_line_stop = rand_x_stop(line_idx); % even layers horizontal

z_start = linspace(perc_val(z, 0.15), perc_val(z, 0.9), n_lowperm_layers);

z_stop = z_start + perc_val(z, 0.02);

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
                        & z == min(z))); % only select ONE horizontal patch
    theta0 = pi/4;
    theta = linspace(theta0, pi - theta0, num_x_cells);
    r = (rand_x_stop(i) - rand_x_start(i))/2;
    z_anticline_start = -r*sin(theta) + z_start(i) + r/2*(1+cos(pi/2-theta0)); % -r*sin to get anticline (since z positive downwards)
    z_anticline_stop = -r*sin(theta) + z_stop(i) + r/2*(1+cos(pi/2-theta0)); % + r*cos(pi/4) to shift cap to bottom of spherical cap
    
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
swr = 0.15;
sor = 0.2;

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
                           'smin', [swr, sor]);
             
s = linspace(0,1,100);                       
krW = fluid.krW(s).';
krO = fluid.krO(1-s).';
invalid_sw = find(s>1-sor);
invalid_so = find(s<swr);
krW(invalid_sw) = nan;
krO(invalid_so) = nan;

plot(s, krW, 'blue', s, krO, 'green');
xlabel('Water saturation');
title('Relative permeability curves');
legend('krW', 'krO', 'Location', 'east');
hold off
                       
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
%rate = 40*meter^3/day; % 25 | 3
bhp = 40*barsa;
% Put well one cell above bottom to avoid it interacting with bottom BC
well_h = 1; % cell perforations in vertical direction %3
perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/10);
W = addWell([], G, rock, perforation_idx, ...
            'Type', 'bhp', 'Val', bhp, ...
            'InnerProduct', 'ip_tpf', ...
            'Radius', 0.1, 'Dir', 'x', ...
            'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
            'Name', 'P1');
 
% W = verticalWell([], G, rock, fix(nx/2), 1, linspace(nz, nz-2, 3), ...
%                  'InnerProduct', 'ip_tpf', ...
%                  'Type', 'bhp', 'Val', bhp, ...
%                  'Radius', 0.1, 'Name', 'I_bot', ...
%                  'Comp_i', [0, 1], 'Sign', 1); % [0 1] composition => co2           
%                      
%% Plot grid
clf;
f1 = figure(1);
plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, perm, 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn); colorbar('southoutside');
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
top_cells = G.cells.indexMap(z < min(z(all_lowperm_cells)));

% hydrostatic in interior right
interior_right_cells = G.cells.indexMap(z >= min(z(all_lowperm_cells)) & x == max(x));
%interior_left_cells = G.cells.indexMap(z >= min(z(all_lowperm_cells)) & x == min(x));
[interior_right_faces, interior_right_fp] = gridCellFaces(G, interior_right_cells);
%[interior_left_faces, interior_left_fp] = gridCellFaces(G, interior_left_cells);
interior_east_edges = interior_right_faces(2:6:numel(interior_right_faces)); % pick every EAST element
%interior_west_edges = interior_left_faces(1:6:numel(interior_left_faces)); % pick every WEST element

pz = fluid.rhoWS * norm(gravity) * uniquetol(z(z >= min(z(all_lowperm_cells)))); % hydrostatic pressure in low-perm region
bc = addBC(bc, interior_east_edges, 'pressure', pz, 'sat', [1 0]);
%bc = addBC(bc, interior_west_edges, 'pressure', pz, 'sat', [1 0]);

dt = rampupTimesteps(12000*day(), 30*day(), 8);
n_steps = numel(dt);

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
% Shut off well halfway
%schedule.step.control(n_steps/2:n_steps) = 2;

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
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f2, strcat(plot_dir, '/sat_0'), 'png');

%% Compute solutions
[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

dt_plot = 10; % plot every x'th iteration
vol_ratios = ones(numel(states)+1, 1);

for i=1:numel(states)
    t = t+dt(i);
    S = states{i}.s(:,1);
    assert(max(S) < 1+eps && min(S) > -eps);
    
    p = reshape(states{i}.pressure, [max(ii), max(kk)]);
    p_mean(:,i) = mean(p, 1).';
    
    figure(2);
    plotCellData(G, S, 'EdgeColor', 'none');
    plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
    colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
    title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow, pause(0.2)  
    
    if ~mod(i, dt_plot)
       filename_f2 = sprintf(strcat(plot_dir, '/sat_%d'), i);
       saveas(f2, filename_f2, 'png');
    end
   
    [tot_vol, top_vol] = Case2Functions.Co2MassRatio(G, dims, gridsize, top_cells, 1-S, rock, fluid);
    vol_ratio = 1 - top_vol/tot_vol; 
    vol_ratios(i+1) = vol_ratio;
end

vol_ratio_filename = sprintf(strcat(data_dir, '/vol_ratio_layers_%d_seed_%d.mat'), n_lowperm_layers, seed.Seed);
save(vol_ratio_filename, 'vol_ratios'); % save to compare for different nr low-perm layers

%% Compare volume ratios for different nr low-perm layers
data_files = dir(strcat(data_dir, '/vol_ratio_layers_*.mat'));

%regexp(data_files(i).name, 'layers_\d+', 'match');

run = {};
used_lab = {};
f3 = figure(3); % only one plot per unique num layer
fi = 3;

for i=1:numel(data_files)
   run{i} = load(strcat(data_files(i).folder, '\', data_files(i).name), 'vol_ratios');
   lab = regexp(data_files(i).name, 'layers_\d+', 'match');
   
   if ~any(ismember(used_lab, lab{1})) % this num layer not plotted yet
      used_lab = cat(2, used_lab, lab{1}); % add to list of used labels
      lab_nr = regexp(lab{1}, '\d+', 'match');
      
      figure(3); % plot only first occurence of this num-layer
      plot(run{i}.vol_ratios, '-o', 'DisplayName', strrep(lab{1}, '_', ' '), 'MarkerSize', 4);
      hold on;
      
      if i ~= 1 % new num layers reached - plot all graphs from previous num layers
          figure(fi)       
          xlabel('Time step');
          drawnow
          hold off
          saveas(fig, vol_layer_filename, 'png');          
      end
      
      fi = fi + 1; % new figure (for new layer number)
      fig = figure(fi);
      figure(fi);
      plot(run{i}.vol_ratios, '-', 'DisplayName', strrep(lab{1}, '_', ' '), 'LineWidth', 1.5);
      vol_layer_filename = sprintf(strcat(plot_base_dir, '/vol_ratio_%s_layers_', date), lab_nr{1});
      title(sprintf('Volume ratio of CO2 trapped in interior - %s layers', lab_nr{1})); 
      hold on;
      
   else
       figure(fi);
       % plot all runs for each unique num layers
       plot(run{i}.vol_ratios, '-', 'DisplayName', strrep(lab{1}, '_', ' '), 'LineWidth', 1.5);
       hold on;
   end
end

% Plot for last num layers must be saved after loop
figure(fi)       
xlabel('Time step');
drawnow
hold off
saveas(fig, vol_layer_filename, 'png'); 

% Save plot for each unique num layer
figure(3);
legend('Location', 'northeast');
xlabel('Time step');
title('Ratio of total CO2 volume residing in low-permeable layers');
drawnow
saveas(f3, strcat(plot_base_dir, '/vol_ratio_', date), 'png');
hold off

%% TESTING
tt = {'hei', 'nh5', 'gt9', 're_4', 'uyt4', 'hei'};
nn = {1, 4, 3, 7, 1, 5};
ss = [];
numel(tt);
for i=1:numel(tt)
   ss = cat(1, ss, [tt(i), nn{i}]);
end