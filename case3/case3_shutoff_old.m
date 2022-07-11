%% 2D CO2 migration - shutoff with leakage
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui
ROOTDIR = strrep(ROOTDIR, '\', '/');
seed = rng();
seed.Seed = 1212; % Must set here, otherwise not updated

%% Define 2D grid
nx = 100; ny = 1; nz = 50; % 100 50
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
n_lowperm_layers = 60;

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case3/plots_shutoff_leakage');
data_dir = strcat(ROOTDIR, '../summer_sintef/case3/data_shutoff_leakage');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_lowperm_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

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

f1 = figure(1);
plot(s, krW, 'blue', s, krO, 'green', 'LineWidth', 1.5);
xlabel('Water saturation');
title('Relative permeability curves');
legend('krW', 'krO', 'Location', 'east');
saveas(f1, strcat(plot_dir, '/relperm'), 'png');
hold off
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, G.cells.num)';
pc_vals = UtilFunctions.Pc(dummy_Sw, 0, 100*barsa, 500*barsa, 2);

region_table = {[dummy_Sw, zeros(numel(dummy_Sw), 1)], [dummy_Sw, pc_vals]}; % container for pc values in each region
region_idx = {setdiff(G.cells.indexMap, all_lowperm_cells), all_lowperm_cells}; % region to interpolate (rest, lowperm)
%fluid_pc.pcOW = @(S, varargin) interpReg(region_table, S, reg_idx); % in both regions, interpolate over saturations S
% or with same function in entire domain
%fluid_pc.pcOW = UtilFunctions.Pc(dummy_Sw, 0, 100*barsa, 500*barsa, 2);
% or using Leverett-J pc to account for differences in permeability

%% Horizontal well
rate = 6*meter^3/day(); % 25 | 3
%bhp = 30*barsa;
% Put well one cell above bottom to avoid it interacting with bottom BC
well_h = 1; % cell perforations in vertical direction
perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/10);
W = addWell([], G, rock, perforation_idx, ...
            'Type', 'rate', 'Val', rate, ...
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
f1 = UtilFunctions.fullsizeFig(1);
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
[interior_right_faces, interior_right_fp] = gridCellFaces(G, interior_right_cells);
interior_east_edges = interior_right_faces(2:6:numel(interior_right_faces)); % pick every EAST element
%interior_west_edges = interior_left_faces(1:6:numel(interior_left_faces)); % pick every WEST element

pz = fluid.rhoWS * norm(gravity) * uniquetol(z(z >= min(z(all_lowperm_cells)))); % hydrostatic pressure in low-perm region
bc = addBC(bc, interior_east_edges, 'pressure', pz, 'sat', [1 0]);

tot_time = 16000*day();
dt = rampupTimesteps(tot_time, 75*day(), 10);
n_steps = numel(dt);

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

%% Shut off well just before halfway
inj_stop = fix(0.3*n_steps);
schedule.control(2) = schedule.control(1); % create second well
schedule.control(2).W.status = 0; % shut off second well
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway

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

f2 = UtilFunctions.fullsizeFig(2); % to hold saturations

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

leaked_ratios = zeros(numel(states)+1, 1);

residual_vol = zeros(numel(states)+1, 1);
residual_ratios = zeros(numel(states)+1, 1);
S_gt_sor = false(numel(states)+1, G.cells.num); % True if oil saturation has reached greater than sor
buff = 0.01; % buffer for residual saturation

for i=1:numel(states)
    t = t+dt(i);
    S = states{i}.s(:,1);
    assert(max(S) < 1+eps && min(S) > -eps);
    
    p = reshape(states{i}.pressure, [max(ii), max(kk)]);
    p_mean(:,i) = mean(p, 1).';
    
    set(f2, 'visible', 'off');
    plotCellData(G, S, 'EdgeColor', 'none');
    plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
    colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
    title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
    axis equal tight;
    view([0, 0]);
    
    
    if ~mod(i, dt_plot)
       set(f2, 'visible', 'on');
       filename_f2 = sprintf(strcat(plot_dir, '/sat_%d'), i);
       saveas(f2, filename_f2, 'png');
    end
   
    [tot_vol, top_vol] = UtilFunctions.Co2VolumeRatio(G, top_cells, 1-S, rock, fluid);
    simulation_vol = min(rate*t, rate*sum(dt(1:inj_stop-1))); % Total CO2 volume capped at point where injection stops   
    vol_ratios(i+1) = (tot_vol - top_vol)/simulation_vol; 
    % Leaked CO2   
    leaked_ratios(i+1) = (simulation_vol - tot_vol) / simulation_vol; % leaked_volume / simulation_vol
    
    % Residually trapped CO2
    S_gt_sor(i+1, :) = ((1-S) > sor);   
    if i >= 10 % avoid plotting outliers
        residual_vol(i+1) = UtilFunctions.Co2Residual(G, 1-S, sor, buff, S_gt_sor, rock);   
        residual_ratios(i+1) = residual_vol(i+1) / simulation_vol; % Ratio of total injected volume trapped in interior domain
    end
end

vol_ratio_filename = sprintf(strcat(data_dir, '/vol_ratio_layers_%d_seed_%d.mat'), n_lowperm_layers, seed.Seed);
leaked_ratio_filename = sprintf(strcat(data_dir, '/leaked_ratio_layers_%d_seed_%d.mat'), n_lowperm_layers, seed.Seed);
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_ratio_layers_%d_seed_%d.mat'), n_lowperm_layers, seed.Seed);
save(vol_ratio_filename, 'vol_ratios'); % save to compare for different nr low-perm layers
save(leaked_ratio_filename, 'leaked_ratios');
save(residual_ratio_filename, 'residual_ratios');

%% Plot interior volume ratio, leakage ratio, residual ratio
vol_files = dir(strcat(data_dir, '/vol_ratio_layers_*.mat'));
leaked_files = dir(strcat(data_dir, '/leaked_ratio_layers_*.mat'));
residual_files = dir(strcat(data_dir, '/residual_ratio_layers_*.mat'));

n_files = numel(vol_files);
used_lab = {};
unique_lab_idx = [];
% Find unique num layers
for i=1:n_files
    lab = regexp(vol_files(i).name, '(?<=layers_)\d+', 'match');
    if ~any(ismember(used_lab, lab{1}))
        used_lab = cat(2, used_lab, lab{1}); % add to list of used labels 
        unique_lab_idx = cat(2, unique_lab_idx, i);
    end
end

unique_lab_num = diff(unique_lab_idx);
if numel(unique_lab_num) > 0
    unique_lab_num = cat(2, unique_lab_num, n_files-unique_lab_idx(numel(unique_lab_idx)));
else
    unique_lab_num = cat(2, unique_lab_num, n_files);
end

run = {{}, {}, {}}; % {volume ratio, leakage, residual}
k = 1; % index to access unique num layers

f3 = figure(3); % only one plot per unique num layer
fi = 2;
fj = 3;

from_i = 1; % index to compute mean/var value from
mean_ratios = {};
var_ratios = {};

for i=1:n_files
   load_vol = load(strcat(vol_files(i).folder, '\', vol_files(i).name), 'vol_ratios');
   run{1}{i} = load_vol.vol_ratios; 
   load_leaked = load(strcat(leaked_files(i).folder, '\', leaked_files(i).name), 'leaked_ratios');
   run{2}{i} = load_leaked.leaked_ratios;
   load_residual = load(strcat(residual_files(i).folder, '\', residual_files(i).name), 'residual_ratios');
   run{3}{i} = load_residual.residual_ratios;
   
   seed_v = regexp(vol_files(i).name, '(?<=seed_)\d+', 'match');
   seed_l = regexp(leaked_files(i).name, '(?<=seed_)\d+', 'match');
   seed_r = regexp(residual_files(i).name, '(?<=seed_)\d+', 'match');
            
   if ismember(i, unique_lab_idx) % new num layer 

      if i ~= 1 % new num layers reached - plot all graphs from previous num layers
          figure(fi)       
          xlabel('Time step');
          xlim([0, n_steps+0.03*n_steps])
          legend('Location', 'northwest');
          drawnow         
          hold off
          saveas(fig, res_layer_filename, 'png');
          
          figure(fj)
          xlabel('Time step');
          xlim([0, n_steps+0.03*n_steps])
          lgd = legend('Location', 'northeast'); lgd.NumColumns = unique_lab_num(k);
          drawnow
          hold off
          saveas(fjg, vl_layer_filename, 'png');         
          
          % Compute mean and variance across layers      
          for j=1:3   
             mean_ratios{j} = mean(cell2mat(run{j}(from_i:i-1)), 2); % average for each time step
             var_ratios{j} = var(cell2mat(run{j}(from_i:i-1)), 0, 2); % i-1 since we don't want to include new member
          end         

          figure(3);
          hold on
          yyaxis left
          plot(mean_ratios{3}, 'Color', 'blue', 'DisplayName', strcat('Mean - layers: ', used_lab{k}), 'LineWidth', 1.5);         
          
          yyaxis right
          plot(var_ratios{3}, 'Color', 'red', 'DisplayName', strcat('Var - layers: ', used_lab{k}));
         
          from_i = i; % start-index for next set of num layers
          k = k + 1;
      end
      
      fi = fi + 2; % new figure (for new layer number)    
      fig = figure(fi);
      figure(fi);
      plot(run{3}{i}, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('Seed %s', seed_r{1}));
      res_layer_filename = sprintf(strcat(plot_base_dir, '/res_ratio_%s_layers_', date), used_lab{k});
      title(sprintf('Ratio of CO2 volume residually trapped - %s layers', used_lab{k}));
      hold on;
      
      fj = fj + 2; % new figure (for new layer number)    
      fjg = UtilFunctions.fullsizeFig(fj);
      figure(fj);
      plot(run{1}{i}, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('VR, seed %s', seed_v{1}));
      hold on;
      plot(run{2}{i}, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('LR, seed %s', seed_l{1}));
      vl_layer_filename = sprintf(strcat(plot_base_dir, '/vl_ratio_%s_layers_', date), used_lab{k});
      title({'Volume ratio of CO2 in interior (VR)', sprintf('and leaked ratio (LR) - %s layers', used_lab{k})});
      hold on;
          
   else
      figure(fi);
      % plot all runs for each unique num layers
      plot(run{3}{i}, '-', 'LineWidth', 1.5,  'DisplayName', sprintf('Seed %s', seed_r{1}));
      hold on;
      
      figure(fj);     
      plot(run{1}{i}, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('VR, seed %s', seed_v{1}));
      hold on;
      plot(run{2}{i}, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('LR, seed %s', seed_l{1}));
      hold on;
      
   end
end
    
% Plot for last num layers must be saved after loop
figure(fi)       
xlabel('Time step');
xlim([0, n_steps+0.03*n_steps])
legend('Location', 'northwest');
drawnow
hold off
saveas(fig, res_layer_filename, 'png'); 

figure(fj)       
xlabel('Time step');
xlim([0, n_steps+0.03*n_steps])
lgd = legend('Location', 'northeast'); lgd.NumColumns = unique_lab_num(numel(unique_lab_num));
drawnow
hold off
saveas(fjg, vl_layer_filename, 'png'); 

% Compute mean and variance for last unique num layer
for j=1:3   
    mean_ratios{j} = mean(cell2mat(run{j}(from_i:n_files)), 2); % average for each time step
    var_ratios{j} = var(cell2mat(run{j}(from_i:n_files)), 0, 2);
end

figure(3);
yyaxis left
plot(mean_ratios{3}, 'Color', 'blue', 'DisplayName', strcat('Mean - layers: ', used_lab{k}), 'LineWidth', 1.5);
ylabel('Mean'); 
hold on;

yyaxis right
plot(var_ratios{3}, 'Color', 'red', 'DisplayName', strcat('Var - layers: ', used_lab{k}));
ylabel('Variance');
hold off;

legend('Location', 'northwest');
xlim([0, n_steps+0.03*n_steps])
xlabel('Time step');
title('Ratio of CO2 volume residually trapped');
drawnow
saveas(f3, strcat(plot_base_dir, '/residual_ratio_', date), 'png');

%% TESTING
