%% Case 4 - gaussian relperm distribution & capillary pressure
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');
seed = rng();
seed.Seed = 8854; % Must set here, otherwise not updated

%% Define 2D grid
nx = 75; ny = 1; nz = 75; % 100 50
lx = 800; ly = 1; lz = 450; % 1000 350
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);
x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);

%% Define rock and fluid objects
lowperm = 10*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [G.cells.num 1]);
poro = 0.3;

%% Directories
n_lowperm_layers = 15;
n_imperm_layers = round(n_lowperm_layers/3);

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/plots_pc_LJ_8');
data_dir = strcat(ROOTDIR, '../summer_sintef/case4/data_pc_LJ_8');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_lowperm_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Set up grid with layers
[x_start, x_stop, z_start, ...
    z_stop, line_idx, anticline_idx] ...
  = GridSetup.DefineLayers(G, n_lowperm_layers, n_imperm_layers);

[perm, all_added_layers, trapped_cells] ...
  = GridSetup.GenerateLayers(line_idx, anticline_idx, ... 
                    n_lowperm_layers, n_imperm_layers, ...
                    G, x_start, x_stop, z_start, z_stop, ...
                    perm, lowperm, baseperm);

%% Compute rock+fluid objects
rock = makeRock(G, perm, poro);
T = computeTrans(G, rock, 'Verbose', true);
swr = 0.15;
sor = 0.2;

% fluid = initSimpleADIFluidJfunc('phases', 'WO', ... % [water, GAS] or OIL?
%                            'mu', [1, 0.05]*centi*poise, ... % viscosity
%                            'n',  [2, 2], ... % relperm powers
%                            'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
%                            'smin', [swr, sor], ...
%                            'rock', rock, ...
%                            'J1', 1.149, ...
%                            'J2', 0.1549, ...
%                            'J3', 0, ...
%                            'k1', 0.994, ...
%                            'k2', 65.0);
                       
fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [1.5, 2], ... % relperm powers
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
dummy_K = linspace(min(perm(added_layers{1})), max(perm(added_layers{1})), G.cells.num)';

p_e = 8*barsa; % 0.5
p_cap = 40*barsa; % 3
median_pc = 8*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 0;

if standard_pc
    plot_pc_lim = p_cap;
    fluid.pcOW = @(S, varargin) runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, G);
else
    plot_pc_lim = 5*median_pc;
    fluid.pcOW = @(S, varargin) UtilFunctions.LeverettJ(S, poro, perm, baseperm, median_pc);
end

%fluid.pcOW = @(S, varargin) runLeverettJ_2D(S, dummy_Sw, poro, perm, dummy_K, baseperm, G);

% fluid.pcOW = @(S, varargin) UtilFunctions.LeverettJnew(S, swr, sor, poro, ...
%                                                      perm, baseperm, 2, 1.149, ...
%                                                      0.1549, 0, 0.994, 65.0);

% figure(1);
% plot(dummy_Sw, pc_vals(:, 200:100:800), 'LineWidth', 1.5);
% xlabel('Water saturation');
% xlim([swr, 1]);
% title('Capillary pressure function');
% saveas(f1, strcat(plot_dir, '/cap_pres'), 'png');
% hold off

%% Horizontal well
rate = 10*meter^3/day(); % 25 | 3
%bhp = 30*barsa;
% Put well one cell above bottom to avoid it interacting with bottom BC
well_h = 1; % cell perforations in vertical direction
perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
W = addWell([], G, rock, perforation_idx, ...
            'Type', 'rate', 'Val', rate, ...
            'InnerProduct', 'ip_tpf', ...
            'Radius', 0.1, 'Dir', 'x', ...
            'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
            'Name', 'P1');
                  
%% Plot grid
clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(perm, milli*darcy); % to better illuminate the gaussian perm in low-perm layers
%perm_dummy(setdiff(1:numel(perm_dummy), all_added_layers)) = nan;
plotCellData(G, log10(perm_dummy), 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn);
colorbarHist(log10(perm_dummy(all_added_layers)), [min(log10(perm_dummy)), max(log10(perm_dummy))], 'South', 51);
title('Log of permeability field');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
drawnow
hold off

saveas(f1, strcat(plot_dir, '/perm_field'), 'png');

%% Set up solver
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
disp(model)

%% Boundary conditions and schedule
bc = []; % no-flux as default
top_cells = G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = G.cells.indexMap(z >= min(z(all_added_layers)));

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, G, 'Right', pz, 'sat', [1 0]);

tot_time = 8000*day();
dt = rampupTimesteps(tot_time, 30*day(), 10);
n_steps = numel(dt);

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

%% Shut off well just before halfway
inj_stop = fix(0.3*n_steps);
schedule.control(2) = schedule.control(1); % create second well
schedule.control(2).W.status = 0; % shut off second well
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at injection stop

%% Initial condition
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(G, 100*barsa, [1,0]);
t = 0;
p_mean = zeros(max(kk), n_steps);

lowperm_faces = zeros([size(all_added_layers), 6]); % 6 faces in 3d cartesian coords
for i=1:size(all_added_layers)
    icell_start = G.cells.facePos(all_added_layers(i));
    icell_stop = G.cells.facePos(all_added_layers(i)+1)-1;
    lowperm_faces(i,:) = G.cells.faces(icell_start:icell_stop, 1);
end

f2 = UtilFunctions.fullsizeFig(2); % to hold saturations

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar; caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f2, strcat(plot_dir, '/sat_0'), 'png');

f3 = UtilFunctions.fullsizeFig(3); % to hold cap pressure

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
%pc_test = fluid.pcOW(state.s(:,1), perm);
%fluid.pcOW(state.s(:,1) == 1) = 0;
plotCellData(G, fluid.pcOW(state.s(:,1)), 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colorbar; caxis([0, plot_pc_lim]);
title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/cap_pres_0'), 'png');

%% Compute solutions
[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

dt_plot = cat(1, repmat(5, [fix(numel(states)/4), 1]), ... % first numel(states)/4 states plotted frequently
                 repmat(10, [numel(states)-fix(numel(states)/4), 1])); % subsequent states plotted less frequently
vol_ratios = ones(numel(states)+1, 1);

leaked_ratios = zeros(numel(states)+1, 1);

residual_vol = zeros(numel(states)+1, 1);
residual_ratios = zeros(numel(states)+1, 1);
simulation_vol = zeros(numel(states)+1, 1);

S_gt_sor = false(numel(states)+1, G.cells.num); % True if oil saturation has reached greater than sor
buff = 0.01; % buffer for residual saturation

for i=1:numel(states)
    t = t+dt(i);
    
    S = states{i}.s(:,1);
    assert(max(S) < 1+eps && min(S) > -eps);
    %fluid.pcOW(S == 1) = 0;
    
    p = reshape(states{i}.pressure, [max(ii), max(kk)]);
    p_mean(:,i) = mean(p, 1).';    
    
    if ~mod(i, dt_plot(i))
       figure(2);
       set(f2, 'visible', 'off');
       plotCellData(G, S, 'EdgeColor', 'none');
       plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
       colormap(flipud(winter)); colorbar; caxis([0 1]);
       title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
       axis equal tight;
       view([0, 0]);   
        
       set(f2, 'visible', 'on');      
       filename_f2 = sprintf(strcat(plot_dir, '/sat_%d'), i);      
       saveas(f2, filename_f2, 'png');
       
       figure(3);
       set(f3, 'visible', 'off');
       plotCellData(G, fluid.pcOW(S), 'EdgeColor', 'none');
       plotGrid(G, W.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
       colorbar; caxis([0, plot_pc_lim]);
       title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
       axis equal tight
       view([0, 0])
    
       set(f3, 'visible', 'on');
       filename_f3 = sprintf(strcat(plot_dir, '/cap_pres_%d'), i);
       saveas(f3, filename_f3, 'png');
    end
   
    [tot_vol, top_vol] = UtilFunctions.Co2VolumeRatio(G, top_cells, 1-S, rock, fluid);
    simulation_vol(i+1) = min(rate*t, rate*sum(dt(1:inj_stop-1))); % Total CO2 volume capped at point where injection stops   
    vol_ratios(i+1) = (tot_vol - top_vol)/simulation_vol(i+1); 
    % Leaked CO2   
    leaked_ratios(i+1) = (simulation_vol(i+1) - tot_vol) / simulation_vol(i+1); % leaked_volume / simulation_vol
    
    % Residually trapped CO2
    % assume imbibition and drainage are mutually exclusive 
    % => imbibition starts when injection stops
    if i >= inj_stop-1
        residual_vol(i+1) = UtilFunctions.Co2ResidualTrapped(G, 1-S, sor, rock);   
        residual_ratios(i+1) = residual_vol(i+1) / simulation_vol(i+1); % Ratio of total injected volume trapped in interior domain
    end       
end

% linear interpolation of residual trapping
dt_samples = [1, inj_stop];
res_samples = residual_vol(dt_samples);
dt_interp = 1:inj_stop;
residual_vol(1:inj_stop) = interp1(dt_samples, res_samples, dt_interp, 'linear');
residual_ratios(1:inj_stop) = residual_vol(1:inj_stop) ./ simulation_vol(1:inj_stop);

% store simulation results
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
    unique_lab_num = cat(2, unique_lab_num, n_files-unique_lab_idx(numel(unique_lab_idx))+1);
else
    unique_lab_num = cat(2, unique_lab_num, n_files);
end

run = {{}, {}, {}}; % {volume ratio, leakage, residual}
k = 1; % index to access unique num layers

f4 = figure(4); % only one plot per unique num layer
fi = 3;
fj = 4;

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

          figure(4);
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

figure(4);
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
saveas(f4, strcat(plot_base_dir, '/residual_ratio_', date), 'png');

%% TESTING

Sw = linspace(0, 1, 10).';
Kw = linspace(10*milli*darcy, 100*milli*darcy, 10).';
[Sw, Kw] = ndgrid(Sw, Kw);

%J_vals = UtilFunctions.LeverettJ(Sw, 0.15, 0.2, 0.5, Kw, 2, 1.149, 0.1549, 0, 0.994, 65.0);

%% Functions
function pc = runStandardPc(S, dummy_S, swr, snr, p_e, p_cap, layers, G)
    pc_vals = UtilFunctions.PcNew(dummy_S, swr, snr, p_e, p_cap, 2);

    region_table = {[dummy_S, zeros(numel(dummy_S), 1)], [dummy_S, pc_vals]}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg(region_table, S, region_idx);
end

function pc = runLeverettJ_2D(S, dummy_S, phi, K, dummy_K, K_base, layers, G)
    [grid_Sw, grid_K] = ndgrid(dummy_S, dummy_K);
    pc_vals = UtilFunctions.LeverettJ(grid_Sw, phi, grid_K, K_base);
    
    region_table = {{grid_Sw, grid_K, zeros(size(grid_Sw))}, ...
                     {grid_Sw, grid_K,  pc_vals}}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg2D(region_table, S, K, region_idx);
end