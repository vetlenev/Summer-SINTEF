%% Case 4 - gaussian relperm distribution & capillary pressure
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 1313;
n_layers = 13;
n_cuts = struct('layers13', 2, 'layers26', 3, 'layers40', 4);

%% Define 2D grid
nx = 60; ny = 1; nz = 40; % 100 50
lx = 800; ly = 1; lz = 400; % 1000 350
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

[ii, jj, kk] = gridLogicalIndices(G);
x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);
dz = mean(diff(z)); % average cell spacing z-dir

%% Define rock and fluid objects
lowperm = 10*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [G.cells.num 1]);
perms = {lowperm, baseperm};
poro = 0.3;

%% Directories
n_lowperm_layers = ceil(n_layers/2); 
n_imperm_layers = floor(n_layers/2);

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case5/plots_hysteresis');
data_base_dir = strcat(ROOTDIR, '../summer_sintef/case5/data_hysteresis');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_layers);
data_dir = sprintf(strcat(data_base_dir, '/layers_%d'), n_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

plot_dir_base = strcat(plot_dir, '/base');
plot_dir_hys = strcat(plot_dir, '/hys');
data_dir_base = strcat(data_dir, '/base');
data_dir_hys = strcat(data_dir, '/hys');
dir_exists = mkdir(plot_dir_base) & mkdir(data_dir_base) & mkdir(data_dir_hys) & mkdir(plot_dir_hys);

this_lp = regexp(data_dir, 'lp_\d+', 'match');
this_layer = regexp(data_dir, 'layers_\d+', 'match');

%% Set seed
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Set up grid with layers
grid = GridSetupParams(G, n_lowperm_layers, n_imperm_layers, perms);

[line_idx, anticline_idx] = grid.DefineLayers();

grid.angle = pi/4;
grid.log_min = 0.1; grid.log_max = 10;
grid.nCuts = n_cuts.(strcat('layers', string(n_layers)));

[added_layers, trapped_cells, ...
    line_idxs, anticline_idxs, ...
    anticline_cells_dummy] = grid.GenerateLayers(line_idx, anticline_idx);

all_lowperm_layers = added_layers{2};
all_lowperm_layers = vertcat(all_lowperm_layers{:});
all_imperm_layers = added_layers{1};
all_imperm_layers = vertcat(all_imperm_layers{:});
all_added_layers = [added_layers{:}];
all_added_layers = vertcat(all_added_layers{:});

grid.corr_len_x = 3;
grid.corr_len_z = 3;

[new_trapped_cells] = grid.SetPermeability(added_layers, trapped_cells, ...
                                        line_idxs, anticline_idxs, ...
                                        anticline_cells_dummy);   

s_trapped_imperm = new_trapped_cells{1};
s_trapped_free = new_trapped_cells{2};

%% Compute rock+fluid objects
rock = makeRock(grid.G, grid.Perm, poro);
T = computeTrans(grid.G, rock, 'Verbose', true);

swr = 0.15;
snr = 0.2;

% fluid for primary drainage
fluid_P = initSimpleADIFluid('phases', 'WO', ...
                           'mu', [1, 0.05]*centi*poise, ...
                           'n',  [2, 2], ...
                           'rho', [1000, 650]*kilogram/meter^3, ... 
                           'smin', [swr, 0]);

% fluid standard                     
fluid_base = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
                           'smin', [swr, snr]);
                       
% fluid hysteresis
fluid_hys = fluid_base;

lambda = 2;
sni_max = 1-swr;
snr_max = snr;

krn_PD =  @(s) fluid_P.krO(s); % primary drainage
krn_PI = @(s) fluid_base.krO(s); % primary imbibition

%krn_PI =  @(s) Hysteresis.Killough(s.val, sni_max, sni_max, snr_max, lambda, krn_PD); % primary imbibition
%fluid_hys.krO = @(s, varargin) Hysteresis.KilloughOld(s, varargin{1}, sni_max, snr_max, lambda, krn_PD); % varargin{1} = sni

fluid_hys.krO = @(s, sMax) Hysteresis.Killough(s, sMax, sni_max, snr_max, krn_PD, krn_PI);

             
%% Compare relperms
s = linspace(0,1,100);                       
krW = fluid_base.krW(s).';
krO = fluid_base.krO(1-s).';
krW(s>1-snr) = nan;
krO(s<swr) = nan;

krW_P = fluid_P.krW(s).';
krO_P = fluid_P.krO(1-s).';
krW_P(s>1) = nan;
krO_P(s<swr) = nan;

sMax = 0.7;
krW_H = fluid_hys.krW(s).';
krO_H = fluid_hys.krO(1-s, sMax);
krW_H(s>1-snr) = nan;
krO_H(s<swr) = nan;

f10 = figure(10);
plot(s, krO_P, 'blue', s, krO, 'red', 'LineWidth', 1.5);
hold on
plot(s, krO_H, '--g', 'LineWidth', 1.5);
xlabel('Water saturation');
title('Relative permeability hysteresis');
legend('krO_{PD}', 'krO_{PI}', 'krO_H', 'Location', 'east');
drawnow
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, grid.G.cells.num)';
p_e = 1*barsa;
p_cap = 5*barsa;

median_pc = 1*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 1;

if standard_pc
    plot_pc_lim = 2*p_e;
    fluid_base.pcOW = @(S, varargin) Capillary.runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, grid.G);
    fluid_hys.pcOW = @(S, varargin) Capillary.runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, grid.G);
else
    plot_pc_lim = 5*median_pc;
    fluid_base.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, perm, baseperm, median_pc);
    fluid_hys.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, grid.Perm, baseperm, median_pc);
end

%% Dummy well
well_h = 1; % cell perforations in vertical direction
perforation_idx = grid.G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
dummy_rock = makeRock(grid.G, grid.Perm, poro);
dummyW = addWell([], grid.G, dummy_rock, perforation_idx, ...
        'Type', 'rate', 'Val', 1*meter^3/day(), ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');
    
%% Plot permeability field
all_trapped_cells = horzcat(trapped_cells{:});

clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(grid.G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(grid.Perm, milli*darcy);
plotCellData(grid.G, log10(perm_dummy), 'EdgeColor', 'none');
plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(jet));
colorbarHist(log10(perm_dummy(all_added_layers)), [min(log10(perm_dummy)), max(log10(perm_dummy))], 'South', 51);
title('Log of permeability field');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
drawnow
hold off

saveas(f1, strcat(plot_dir, '/perm'), 'png');        

%% Plot structurally trapped cells
f2 = UtilFunctions.fullsizeFig(2);
plotGrid(grid.G, all_lowperm_layers, 'FaceColor', 'yellow', 'EdgeColor', 'black', 'EdgeAlpha', 0.2, 'DisplayName', 'Lowperm layers');
plotGrid(grid.G, all_imperm_layers, 'FaceColor', 'red', 'EdgeColor', 'black', 'EdgeAlpha', 0.2, 'DisplayName', 'Imperm layers');
plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none', 'DisplayName', 'Well');
plotGrid(grid.G, s_trapped_free, 'FaceColor', [0.25, 0.5, 0.25], 'EdgeColor', 'none', 'DisplayName', 'Lowperm trapped');
plotGrid(grid.G, s_trapped_imperm, 'FaceColor', [0.5, 0, 0.5], 'EdgeColor', 'none', 'DisplayName', 'Imperm trapped');
title('Structurally trapped CO2');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
[hleg, hobj] = legend('Location', 'southoutside', 'Orientation', 'horizontal');
textobj = findobj(hobj, 'type', 'text');
set(hleg, 'position', [0.1, 0.06, 0.8, 0.05]);
set(textobj, 'fontsize', 12);
legend('Location', 'southoutside', 'Orientation', 'horizontal');
drawnow
hold off

saveas(f2, strcat(plot_dir, '/struct_trapped'), 'png');  


%% Set up solver
gravity reset on
fluids = {fluid_base, fluid_hys};
model_base = TwoPhaseOilWaterModel(grid.G, rock, fluid_base);
model_hys = TwoPhaseOilWaterModel(grid.G, rock, fluid_hys);
disp(model_base)

%% Tailor the hysteresis model
model_hys = model_hys.validateModel;

model_hys.FlowPropertyFunctions.RelativePermeability = ...
                            RelativePermeabilityHys(model_hys);
                        
disp(model_hys)               

%% Boundary conditions, schedule and initial conditions
top_cells = grid.G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = grid.G.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid_base.rhoWS * norm(gravity) * min(z);
bc = pside(bc, grid.G, 'Top', p_top, 'sat', [1 0]);

pz = fluid_base.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, grid.G, 'Right', pz, 'sat', [1 0]);

tot_time = 1000*365*day();
%dt = rampupTimesteps(tot_time, 1000*day(), 10);
dt = rampupTimesteps(round(tot_time/10), 250*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/10), 1000*day(), 0)']';
num_dt = numel(dt);

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

inj_stop_rate = 0.05;
optimal_rates = struct('layers13', 0.25, 'layers26', 0.35, 'layers40', 0.45); % corresponding to approx 5% leakage
rate = optimal_rates.(strcat('layers',string(n_layers)))*meter^3/day();

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(grid.G, 100*barsa, [1,0]);
t = 0;
fig_i = 3;
fig_sat = {};
fig_cap = {};

for i=1:numel(fluids)
    fig_sat{i} = UtilFunctions.fullsizeFig(fig_i); % to hold saturations
    plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
    plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
    plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
    colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
    title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow
    
    fig_cap{i} = UtilFunctions.fullsizeFig(fig_i+1); % to hold cap pressure    
    plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
    plotCellData(G, fluids{i}.pcOW(state.s(:,1)), 'EdgeColor', 'none');
    plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
    colorbar; caxis([0, plot_pc_lim]);
    title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow
    
    fig_i = fig_i + 2;
end

%% Run experiment
[max_volume_hys, leaked_boundary_hys, ...
        structural_utilized_hys, ...
        states_hys, categorized_vols_hys] = RunSimComputeTrapping(grid, rock, rate, state, model_hys, ...
                                                          s_trapped_imperm, s_trapped_free, ...
                                                          swr, snr, dt, bc, inj_stop_rate);

[max_volume_base, leaked_boundary_base, ...
        structural_utilized_base, ...
        states_base, categorized_vols_base] = RunSimComputeTrapping(grid, rock, rate, state, model_base, ...
                                                          s_trapped_imperm, s_trapped_free, ...
                                                          swr, snr, dt, bc, inj_stop_rate);                                                                                                           

fluid_states = {states_base, states_hys};
                                                                     
%% Show optimal solution found                                                       
dt_plot = cat(1, repmat(10, [fix(numel(dt)/10), 1]), ...
                 repmat(30, [numel(dt)-fix(numel(dt)/10), 1]));

sim_type = {'base', 'hys'}; % base, hysteresis
fig_i = 3;
             
for k=1:numel(fluids)
    
    states = fluid_states{k};
    fluid = fluids{k};
    t = 0;    
    
    for i=1:numel(states)
        t = t + dt(i);   

        if ~mod(i, dt_plot(i))
            figure(fig_i)
            set(fig_sat{k}, 'visible', 'off');
            plotCellData(grid.G, states{i}.s(:,1), 'EdgeColor', 'none');
            plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
            colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
            title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
            axis equal tight;
            view([0, 0]);   

            set(fig_sat{k}, 'visible', 'on');
            filename_f3 = sprintf(strcat(plot_dir, '/%s/sat%d_%d'), sim_type{k}, k, i);
            saveas(fig_sat{k}, filename_f3, 'png');

            figure(fig_i+1);
            set(fig_cap{k}, 'visible', 'off');
            plotCellData(grid.G, fluid_base.pcOW(states{i}.s(:,1)), 'EdgeColor', 'none');
            plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
            colorbar; caxis([0, plot_pc_lim]);
            title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
            axis equal tight
            view([0, 0])

            set(fig_cap{k}, 'visible', 'on');
            filename_f4 = sprintf(strcat(plot_dir, '/%s/cap_pres%d_%d'), sim_type{k}, k, i);
            saveas(fig_cap{k}, filename_f4, 'png');
        end 
    end
    
    fig_i = fig_i + 2;
end


%% Store volume categories
categorized_vols = {categorized_vols_base, categorized_vols_hys};
structural_utilized = {structural_utilized_base, structural_utilized_hys};

for k=1:numel(sim_type)

    residual_ratio_filename = sprintf(strcat(data_dir, '/%s/residual_seed_%d.mat'), sim_type{k}, seed.Seed);
    struct_ratio_filename = sprintf(strcat(data_dir, '/%s/struct_seed_%d.mat'), sim_type{k}, seed.Seed);
    free_ratio_filename = sprintf(strcat(data_dir, '/%s/free_seed_%d.mat'), sim_type{k}, seed.Seed);
    exit_ratio_filename = sprintf(strcat(data_dir, '/%s/exit_seed_%d.mat'), sim_type{k}, seed.Seed);
    struct_util_filename = sprintf(strcat(data_dir, '/%s/util_seed_%d.mat'), sim_type{k}, seed.Seed);

    residual_vol = categorized_vols{k}{1};
    structural_vol = categorized_vols{k}{2};
    free_vol = categorized_vols{k}{3};
    exit_vol = categorized_vols{k}{4};
    struct_utilized = structural_utilized{k};

    save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
    save(struct_ratio_filename, 'structural_vol');
    save(free_ratio_filename, 'free_vol');
    save(exit_ratio_filename, 'exit_vol');
    save(struct_util_filename, 'struct_utilized');
end

%% Read data and store in structs
structural_struct = struct; % struct to hold structurally trapped values for each leaked percentage and num layer
residual_struct = struct;
free_struct = struct;
exit_struct = struct;
util_struct = struct;

layer_folders = dir(strcat(data_base_dir, '/layers_*'));
layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
this_layer = strcat('layers_', string(n_layers));

num_layers = numel(layer_folders);

used_lab = {};

for j=1:num_layers
    
    lab_layers = regexp(layer_folders(j).name, 'layers_\d+', 'match');      
    used_lab = cat(1, used_lab, lab_layers{1});
    
    for i=1:numel(sim_type)     

        structural_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/', sim_type{i}, '/struct_*.mat'));
        free_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/', sim_type{i}, '/free_*.mat'));
        residual_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/', sim_type{i}, '/residual_*.mat'));
        exit_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/', sim_type{i}, '/exit_*.mat'));

        util_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/', sim_type{i}, '/util_*.mat'));

        for k=1:numel(free_files) % every seed
            get_seed = regexp(free_files(k).name, '\d+', 'match');
            get_seed = strcat('s', get_seed{1});

            load_structural = load(strcat(structural_files(k).folder, '\', structural_files(k).name), 'structural_vol');
            load_structural = load_structural.structural_vol;
            load_free = load(strcat(free_files(k).folder, '\', free_files(k).name), 'free_vol');
            load_free = load_free.free_vol;
            load_residual = load(strcat(residual_files(k).folder, '\', residual_files(k).name), 'residual_vol');
            load_residual = load_residual.residual_vol;
            load_exit = load(strcat(exit_files(k).folder, '\', exit_files(k).name), 'exit_vol');
            load_exit = load_exit.exit_vol;

            load_util = load(strcat(util_files(k).folder, '\', util_files(k).name), 'struct_utilized');
            load_util = load_util.struct_utilized;
            
            structural_struct.imperm.(lab_layers{1}).(sim_type{i}).(get_seed) = load_structural.imperm;             
            structural_struct.lowperm.(lab_layers{1}).(sim_type{i}).(get_seed) = load_structural.lowperm;
                                                                                                  
            free_struct.(lab_layers{1}).(sim_type{i}).(get_seed) = load_free;   
            residual_struct.(lab_layers{1}).(sim_type{i}).(get_seed) = load_residual;
            exit_struct.(lab_layers{1}).(sim_type{i}).(get_seed) = load_exit;
          
            util_struct.(sim_type{i}).(lab_layers{1}).(get_seed) = load_util;

        end
    end
end


%% Copmarison of trapping inventory, for each num layer

mean_trapping = struct;

log_time = vertcat(0, cumsum(dt)*second()/year());

for i=1:num_layers
    layer = used_lab{i};
    layer_int = regexp(layer, '\d+', 'match');
    layer_int = str2double(layer_int{1});
    
    all_mean_traps = {};
    
    for j=1:numel(sim_type)       
       mean_trapping.(layer).(sim_type{j}).structural.imperm = mean(struct2array(structural_struct.imperm.(layer).(sim_type{j})), 2);
       mean_trapping.(layer).(sim_type{j}).structural.lowperm = mean(struct2array(structural_struct.lowperm.(layer).(sim_type{j})), 2);
       mean_trapping.(layer).(sim_type{j}).residual = mean(struct2array(residual_struct.(layer).(sim_type{j})), 2);
       mean_trapping.(layer).(sim_type{j}).free = mean(struct2array(free_struct.(layer).(sim_type{j})), 2);
       mean_trapping.(layer).(sim_type{j}).exit = mean(struct2array(exit_struct.(layer).(sim_type{j})), 2);
       
       stacked_trapping_mean = [mean_trapping.(layer).(sim_type{j}).structural.imperm.'; ...
                                mean_trapping.(layer).(sim_type{j}).structural.lowperm.'; ...
                                mean_trapping.(layer).(sim_type{j}).residual.'; ...
                                mean_trapping.(layer).(sim_type{j}).free.'; ...
                                mean_trapping.(layer).(sim_type{j}).exit.'].';
                            
        all_mean_traps = cat(1, all_mean_traps, stacked_trapping_mean);
                           
        fig_mean = figure(fig_i);
        figure(fig_i)
        area(log_time, stacked_trapping_mean)
        set(gca, 'xscale', 'log')

        xlim([1, log_time(end)]) % plot from 1 year to end
        xlabel('Time step')
        ylabel('Volume (m^3)')
        title({sprintf('Model: %s', sim_type{j}), sprintf('Mean trapping distribution: %d layers', layer_int)})  
        legend('Structural permanent', 'Structural temporary', ...
                'Residual', 'Free plume', 'Exited', 'Location', 'west')  
        saveas(fig_mean, sprintf(strcat(plot_base_dir, '/%s/%s/mean_trapping'), layer, sim_type{j}), 'png'); 
           
        fig_i = fig_i + 1;
    end
    
    fig_diff = figure(fig_i);
    figure(fig_i)
    
    mean_diff = abs(minus(all_mean_traps{:}));
    
    area(log_time, mean_diff)
    set(gca, 'xscale', 'log')

    xlim([1, log_time(end)]) % plot from 1 year to end
    xlabel('Time step')
    ylabel('Volume (m^3)')
    title({'Hysteresis vs. no hysteresis.', sprintf('Abs difference in trapping inventory (%d layers)', layer_int)})  
    legend('Structural permanent', 'Structural temporary', ...
            'Residual', 'Free plume', 'Exited', 'Location', 'west')  
    saveas(fig_diff, sprintf(strcat(plot_base_dir, '/%s/diff_trapping'), layer), 'png'); 

    fig_i = fig_i + 1;

end

%% Plot mean+var utilized struct traps
clr = {'blue', 'red', 'green', 'magenta', 'orange'};
linSmean = {'-^', '--o', '--x', '--p'};
linWidth = [5, 3, 1.5, 1.0]; 
markS = [6, 6, 10, 5];

fig_util = figure(fig_i);
figure(fig_i)

for i=1:numel(sim_type)  
    mean_util = [];
    std_util = [];
    
    for j=1:num_layers
        layer = used_lab{j};        

        mean_util = cat(1, mean_util, mean(struct2array(util_struct.(sim_type{i}).(layer)), 2));
        std_util = cat(1, std_util, std(struct2array(util_struct.(sim_type{i}).(layer)), 0, 2));   
    end
    
    errorbar(1:num_layers, mean_util, std_util, linSmean{i}, ...
                'CapSize', 18, 'Color', clr{i}, 'MarkerFaceColor', clr{i}, 'LineWidth', 1.5, ...
                'MarkerSize', markS(i), 'DisplayName', sim_type{i});
    hold on
end

xticks(1:num_layers);
xticklabels(strrep(used_lab, 'layers_', ''));
xlabel('Number of low+imperm layers')
ylabel('Utilized (%)')
title({'Percentage of permanent structural traps utilized'})
legend()
drawnow

saveas(fig_util, strcat(plot_base_dir, '/struct_utilized'), 'png');

%% TESTING 2


