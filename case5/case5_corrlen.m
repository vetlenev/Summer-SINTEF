%% Case 5 - dependency of trapping on various parameters (here: corr_len)
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 7755;

%% Set simulation settings
param_name = 'corrlen';
n_layers = 13;

%% Define 2D grid
nx = 60; ny = 1; nz = 40; % 100 50
lx = 800; ly = 1; lz = 400; % 1000 350
dims = [nx ny nz];
gridsize = [lx, ly, lz]*meter;
global G; % global to be accessed inside functions
G = cartGrid(dims, gridsize);
G = computeGeometry(G);

x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);

%% Define rock and fluid objects
lowperm = 1*milli*darcy;
baseperm = 100*milli*darcy;
perms = {lowperm, baseperm};

poro = 0.3;
gravity reset on

%% Set parameters
%corr_len_x = mean(grid.xStop.lowperm - grid.xStart.lowperm) ./ params;
corr_len_x = [27, 11, 5, 3, 1];
%corr_len_z = mean(grid.zStop.lowperm - grid.zStart.lowperm) ./ (params*0.5);
corr_len_z = [19, 7, 5, 3, 1];

params = 1:numel(corr_len_x);
plot_sat = params(1:2:end);

n_cuts = struct('layers13', 2, 'layers26', 3, 'layers40', 4);

%% Directories
n_lowperm_layers = ceil(n_layers/2); 
n_imperm_layers = floor(n_layers/2);

plot_base_dir = strcat(ROOTDIR, sprintf('../summer_sintef/case5/plots_%s', param_name));
data_base_dir = strcat(ROOTDIR, sprintf('../summer_sintef/case5/data_%s', param_name));

plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d'), n_layers);
data_dir = sprintf(strcat(data_base_dir, '/layers_%d'), n_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Set seed
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

%% Set up grid object with layers
grid = GridSetupParams(G, n_lowperm_layers, n_imperm_layers, perms);

[line_idx, anticline_idx] = grid.DefineLayers();

grid.angle = pi/4;
grid.log_min = 0.05; grid.log_max = 100;
grid.nCuts = n_cuts.(strcat('layers', string(n_layers)));
% grid.corr_len_x = mean(x_stop.lowperm - x_start.lowperm) / 100;
% grid.corr_len_z = mean(z_stop.lowperm - z_start.lowperm) / 10;

[added_layers, trapped_cells, ...
    line_idxs, anticline_idxs, ...
    anticline_cells_dummy] = grid.GenerateLayers(line_idx, anticline_idx);

all_lowperm_layers = added_layers{2};
all_lowperm_layers = vertcat(all_lowperm_layers{:});
all_imperm_layers = added_layers{1};
all_imperm_layers = vertcat(all_imperm_layers{:});
all_added_layers = [added_layers{:}];
all_added_layers = vertcat(all_added_layers{:});
                                                    
%% Compute fluid objects
swr = 0.15;
snr = 0.2;

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.05]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 650]*kilogram/meter^3, ... % densities: [water, CO2]
                           'smin', [swr, snr]);
             
s = linspace(0,1,100);                       
krW = fluid.krW(s).';
krO = fluid.krO(1-s).';
invalid_sw = find(s>1-snr);
invalid_so = find(s<swr);
krW(invalid_sw) = nan;
krO(invalid_so) = nan;

f1 = figure(1);
plot(s, krW, 'blue', s, krO, 'green', 'LineWidth', 1.5);
xlabel('Water saturation');
title('Relative permeability curves');
legend('krW', 'krO', 'Location', 'east');
saveas(f1, strcat(plot_dir, '/relperm'), 'png');
clf;
                       
%% Define capillary pressure
dummy_Sw = linspace(0, 1, grid.G.cells.num)';
p_e = 1*barsa;
p_cap = 4*barsa;

median_pc = 0.5*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 0;

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
    
%% Boundary conditions, schedule
top_cells = grid.G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = grid.G.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, grid.G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, grid.G, 'Right', pz, 'sat', [1 0]);

tot_time = 1000*365*day();
dt = rampupTimesteps(round(tot_time/10), 500*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/10), 1000*day(), 0)']';

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

%% Initial state - same for all parameter values
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(grid.G, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotGrid(grid.G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(grid.G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

%% Simulation params
inj_stop_rate = 0.05;
% start with injection rate corresponding to 1/pv'th of total pore volume
optimal_rates = struct('layers13', 0.25, 'layers26', 0.35, 'layers40', 0.45); % corresponding to approx 5% leakage
rate = optimal_rates.(strcat('layers',string(n_layers)))*meter^3/day();

%% Run simulations
allResults = struct;
params_str = cell(1, numel(params));

for p=1:numel(params)  
    grid.corr_len_x =  corr_len_x(p);
    grid.corr_len_z = corr_len_z(p);       

    [new_trapped_cells] = grid.SetPermeability(added_layers, trapped_cells, ...
                                        line_idxs, anticline_idxs, ...
                                        anticline_cells_dummy);   

    s_trapped_imperm = new_trapped_cells{1};
    s_trapped_free = new_trapped_cells{2};
    
    % Compute rock properties
    rock = makeRock(grid.G, grid.Perm, poro);
    T = computeTrans(grid.G, rock, 'Verbose', true);
    
    % Domain-specific capillary pressure
    if standard_pc
        plot_pc_lim = p_cap;
        fluid.pcOW = @(S, varargin) runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, grid.G);
    else
        plot_pc_lim = 5*median_pc;
        fluid.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, grid.Perm, baseperm, median_pc);
    end
    
    % Define model
    model = TwoPhaseOilWaterModel(grid.G, rock, fluid);       

    % Plot permeability field
    figure(1)
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

    saveas(f1, sprintf(strcat(plot_dir, '/perm_p%d'), p), 'png');        

    % Plot structurally trapped cells
    figure(2)
    clf;
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

    % Run simulations
    [max_volume, leaked_boundary, ...
        structural_utilized, ...
        states, categorized_vols] = RunSimComputeTrapping(grid, rock, rate, state, model, ...
                                                          s_trapped_imperm, s_trapped_free, ...
                                                          swr, snr, dt, bc, inj_stop_rate);
    
    params_str{p} = strcat('p', string(p));
    
    allResults.max_volumes.(params_str{p}) = max_volume;
    allResults.leaked_boundary.(params_str{p}) = leaked_boundary;
    allResults.states.(params_str{p}) = states;
    allResults.categorized_vols.(params_str{p}) = categorized_vols;
    allResults.struct_util.(params_str{p}) = structural_utilized;
end

%% Show solution for one of params
%show_param = params(randi([1 numel(params)]));
for p=plot_sat

    show_param_idx = strcat('p', string(find(params == p)));
    show_states = allResults.states.(show_param_idx);
    show_leaked_boundary = allResults.leaked_boundary.(show_param_idx);

    dt_plot = cat(1, repmat(10, [fix(numel(show_states)/5), 1]), ...
                     repmat(25, [numel(show_states)-fix(numel(show_states)/5), 1]));
    
    t = 0;
    for i=1:numel(states)
        t = t + dt(i);   

        if ~mod(i, dt_plot(i))
            figure(3)
            set(f3, 'visible', 'off');
            plotCellData(grid.G, show_states{i}.s(:,1), 'EdgeColor', 'none');
            plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
            colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
            title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)] ...
                    sprintf('Param: %s', show_param_idx)});
            axis equal tight;
            view([0, 0]);   

            set(f3, 'visible', 'on');
            filename_f3 = sprintf(strcat(plot_dir, '/%s_sat_%d'), show_param_idx, i);
            saveas(f3, filename_f3, 'png');

    %         figure(4);
    %         set(f4, 'visible', 'off');
    %         plotCellData(grid.G, fluid.pcOW(show_states{i}.s(:,1)), 'EdgeColor', 'none');
    %         plotGrid(grid.G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
    %         colorbar; caxis([0, plot_pc_lim]);
    %         title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)], ...
    %                 ['Param value: ', show_param]});
    %         axis equal tight
    %         view([0, 0])
    %        
    %         set(f4, 'visible', 'on');
    %         filename_f4 = sprintf(strcat(plot_dir, '/cap_pres_%d_param_%d'), i, show_param);
    %         saveas(f4, filename_f4, 'png');
        end 
    end

    fprintf('\nLeaked boundary for: %s (x) = %.2f, %s (z) = %.2f\n', param_name, corr_len_x(p), param_name, corr_len_z(p))
    disp(show_leaked_boundary)
end

%% Store volume categories
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_%s_seed_%d.mat'), param_name, seed.Seed);
struct_ratio_filename = sprintf(strcat(data_dir, '/struct_%s_seed_%d.mat'), param_name, seed.Seed);
free_ratio_filename = sprintf(strcat(data_dir, '/free_%s_seed_%d.mat'), param_name, seed.Seed);
exit_ratio_filename = sprintf(strcat(data_dir, '/exit_%s_seed_%d.mat'), param_name, seed.Seed);
%inj_rate_filename = sprintf(strcat(data_dir, '/inj_rate_%s_seed_%d.mat'), param_name, seed.Seed);
struct_util_filename = sprintf(strcat(data_dir, '/util_%s_seed_%d.mat'), param_name, seed.Seed);

residual_vol = struct; structural_vol = struct; 
free_vol = struct; exit_vol = struct;
struct_utilized = struct;

for p=1:numel(params)
    p_idx = strcat('p', string(find(params == params(p))));
    residual_vol.(p_idx) = allResults.categorized_vols.(params_str{p}){1};
    structural_vol.(p_idx) = allResults.categorized_vols.(params_str{p}){2};
    free_vol.(p_idx) = allResults.categorized_vols.(params_str{p}){3};
    exit_vol.(p_idx) = allResults.categorized_vols.(params_str{p}){4};
    %inj_rate = rates(end);
    struct_utilized.(p_idx) = allResults.struct_util.(params_str{p});
end

save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
save(struct_ratio_filename, 'structural_vol');
save(free_ratio_filename, 'free_vol');
save(exit_ratio_filename, 'exit_vol');
%save(inj_rate_filename, 'inj_rate');
save(struct_util_filename, 'struct_utilized');

%% Read data and store in structs
layer_folders = dir(strcat(data_base_dir, '/layers_*'));
layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
this_layer = strcat('layers_', string(n_layers));

num_layers = numel(layer_folders);

structural_struct = struct; % struct to hold structurally trapped values for each leaked percentage and num layer
residual_struct = struct;
free_struct = struct;
exit_struct = struct;
%inj_struct = struct;  
util_struct = struct;

used_lab = {};

for j=1:num_layers
    lab_layers = regexp(layer_folders(j).name, 'layers_\d+', 'match');      
    used_lab = cat(1, used_lab, lab_layers{1});

    structural_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/struct_*.mat'));
    free_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/free_*.mat'));
    residual_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/residual_*.mat'));
    exit_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/exit_*.mat'));
    
    util_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/util_*.mat'));

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

        for p=1:numel(params)
            p_idx = strcat('p', string(find(params == params(p))));
            
            structural_struct.imperm.(lab_layers{1}).(p_idx).(get_seed) = load_structural.(p_idx).imperm;             
            structural_struct.lowperm.(lab_layers{1}).(p_idx).(get_seed) = load_structural.(p_idx).lowperm;
                                                                                                  
            free_struct.(lab_layers{1}).(p_idx).(get_seed) = load_free.(p_idx);   
            residual_struct.(lab_layers{1}).(p_idx).(get_seed) = load_residual.(p_idx);
            exit_struct.(lab_layers{1}).(p_idx).(get_seed) = load_exit.(p_idx);
          
            util_struct.(lab_layers{1}).(p_idx).(get_seed) = load_util.(p_idx);
        end
    end
end


%% Plot mean trapping distributions
fm = 7;
linS = {'-', '--', ':', '-o', '--o'};
mycolors = [0 0.5 1; 0.5 0 1; 0 0.5 0; 1 0.647 0; 1 0 0];
%mycolors = {'blue', 'cyan', 'green', 'yellow', 'red'};             

mean_params = struct;
std_params = struct; % 3 figures (one per layer)
std_layers = struct; % 5 figures (one per param)
      
log_time = vertcat(0, cumsum(dt)*second()/year());

for l=1:num_layers
    layer = used_lab(l);
    layer = layer{1};
    for p=1:numel(params)
        p_idx = strcat('p', string(find(params == params(p))));

        mean_params.(layer).structural.imperm.(p_idx) = mean(struct2array(structural_struct.imperm.(layer).(p_idx)), 2);        
        mean_params.(layer).structural.lowperm.(p_idx) = mean(struct2array(structural_struct.lowperm.(layer).(p_idx)), 2);
        mean_params.(layer).free.(p_idx) = mean(struct2array(free_struct.(layer).(p_idx)), 2);       
        mean_params.(layer).residual.(p_idx) = mean(struct2array(residual_struct.(layer).(p_idx)), 2);       
        mean_params.(layer).exit.(p_idx) = mean(struct2array(exit_struct.(layer).(p_idx)), 2);
       
        stacked_trapping_mean = [mean_params.(layer).structural.imperm.(p_idx).'; ...
                            mean_params.(layer).structural.lowperm.(p_idx).'; ...
                            mean_params.(layer).residual.(p_idx).'; ...
                            mean_params.(layer).free.(p_idx).'; ...
                            mean_params.(layer).exit.(p_idx).'].';           
        
        if strcmp(layer, this_layer) % plot trapping distribution for this layer           
            fig_mean = figure(fm);
            figure(fm)
            %set(fm, 'visible', 'on')
            area(log_time, stacked_trapping_mean)
            set(gca, 'xscale', 'log')

            xlim([1, log_time(end)]) % plot from 1 year to end
            xlabel('Time step')
            ylabel('Volume (m^3)')
            title({sprintf('Mean trapping distribution for %d layers', n_layers), ['X: ', param_name, sprintf(' = %.2f', corr_len_x(p))], ...
                                                                                  ['Z: ', param_name, sprintf(' = %.2f', corr_len_z(p))]})  
            legend('Structural permanent', 'Structural temporary', ...
                    'Residual', 'Free plume', 'Exited', 'Location', 'west')  
            saveas(fig_mean, sprintf(strcat(plot_base_dir, '/layers_%d/mean_trapping_%s'), n_layers, p_idx), 'png'); 

            clf; 
        
            this_seed = strcat('s', string(my_seed));
            stacked_trapping = [structural_struct.imperm.(layer).(p_idx).(this_seed).'; ...
                            structural_struct.lowperm.(layer).(p_idx).(this_seed).'; ...
                            residual_struct.(layer).(p_idx).(this_seed).'; ...
                            free_struct.(layer).(p_idx).(this_seed).'; ...
                            exit_struct.(layer).(p_idx).(this_seed).'].';           
                        
            area(log_time, stacked_trapping)
            set(gca, 'xscale', 'log')

            xlim([1, log_time(end)]) % plot from 1 year to end
            xlabel('Time step')
            ylabel('Volume (m^3)')
            title({sprintf('Trapping distribution for %d layers for seed %d', n_layers, my_seed), ['X: ', param_name, sprintf(' = %.2f', corr_len_x(p))], ...
                                                                                  ['Z: ', param_name, sprintf(' = %.2f', corr_len_z(p))]})  
            legend('Structural permanent', 'Structural temporary', ...
                    'Residual', 'Free plume', 'Exited', 'Location', 'west')       
            fm = fm + 1;                        
        end
              
        %fs = fs + 2;
    end
end

%% Plot standard deviation of trapping distribution
% Plot std of trapping mechanisms across params (for each layer)
fs = fm+1;

for l=1:num_layers
    fig_std = figure(fs);
    figure(fs)
    
    layer = used_lab(l);
    layer = layer{1};
    layers_int = regexp(layer, '\d+', 'match');
    layers_int = str2double(layers_int{1});
    
    std_params.(layer).structural.imperm = std(struct2array(mean_params.(layer).structural.imperm), 0, 2);
    std_params.(layer).structural.lowperm = std(struct2array(mean_params.(layer).structural.lowperm), 0, 2);
    std_params.(layer).residual = std(struct2array(mean_params.(layer).residual), 0, 2);
    std_params.(layer).free = std(struct2array(mean_params.(layer).free), 0, 2);
    std_params.(layer).exit = std(struct2array(mean_params.(layer).exit), 0, 2);
    
    stacked_trapping_std = [std_params.(layer).structural.imperm.'; ...
                        std_params.(layer).structural.lowperm.'; ...
                        std_params.(layer).residual.'; ...
                        std_params.(layer).free.'; ...
                        std_params.(layer).exit.'].';
                    
    area(log_time, stacked_trapping_std)
    set(gca, 'xscale', 'log')

    xlim([1, log_time(end)]) % plot from 1 year to end
    xlabel('Time step')
    ylabel('Volume (m^3)')
    title({sprintf('Standard deviation of trapping distribution for %d layers', layers_int)})  
    legend('Structural permanent', 'Structural temporary', ...
            'Residual', 'Free plume', 'Exited', 'Location', 'west')    
    saveas(fig_std, sprintf(strcat(plot_base_dir, '/layers_%d/std_trapping_all_params'), layers_int), 'png');
    
    fs = fs + 1;
    
    % plot std for each param for THIS num layer
    if strcmp(layer, this_layer)
        for p=1:numel(params)
            fig_std = figure(fs);
            figure(fs)
            p_idx = strcat('p', string(find(params == params(p))));
            
            std_layer.structural.imperm = std(struct2array(structural_struct.imperm.(layer).(p_idx)), 0, 2);
            std_layer.structural.lowperm = std(struct2array(structural_struct.lowperm.(layer).(p_idx)), 0, 2);
            std_layer.residual = std(struct2array(residual_struct.(layer).(p_idx)), 0, 2);
            std_layer.free = std(struct2array(free_struct.(layer).(p_idx)), 0, 2);
            std_layer.exit = std(struct2array(exit_struct.(layer).(p_idx)), 0, 2);

            stacked_trapping_std = [std_layer.structural.imperm.'; ...
                                    std_layer.structural.lowperm.'; ...
                                    std_layer.residual.'; ...
                                    std_layer.free.'; ...
                                    std_layer.exit.'].';

            area(log_time, stacked_trapping_std)
            set(gca, 'xscale', 'log')

            xlim([1, log_time(end)]) % plot from 1 year to end
            xlabel('Time step')
            ylabel('Volume (m^3)')
            title({sprintf('Standard deviation of trapping distribution for %d layers', n_layers), ['X: ', param_name, sprintf(' = %.2f', corr_len_x(p))], ...
                                                                                  ['Z: ', param_name, sprintf(' = %.2f', corr_len_z(p))]})  
            legend('Structural permanent', 'Structural temporary', ...
                    'Residual', 'Free plume', 'Exited', 'Location', 'west')  
            saveas(fig_std, sprintf(strcat(plot_base_dir, '/layers_%d/std_trapping_%s'), n_layers, p_idx), 'png'); 
            
            fs = fs + 1;
        end
    end
end

%% Plot mean+var utilized struct traps
f_util = fs+1;
clr = {'blue', 'red', 'green', 'magenta', 'orange'};
linSmean = {'--^', '--o', '--x', '--p'};
linWidth = [5, 3, 1.5, 1.0]; 
markS = [6, 6, 10, 5];

fig_util = figure(f_util);
figure(f_util)
ticks = {};

for i=1:num_layers
    layer = used_lab{i};
    mean_util = [];
    std_util = [];
    
    for p=1:numel(params)
        p_idx = strcat('p', string(find(params == params(p))));
        if i == 1
            ticks = cat(1, ticks, p_idx);
        end

        mean_util = cat(1, mean_util, mean(struct2array(util_struct.(layer).(p_idx)), 2));
        std_util = cat(1, std_util, std(struct2array(util_struct.(layer).(p_idx)), 0, 2));   
    end
    
    errorbar(1:numel(params), mean_util, std_util, linSmean{i}, ...
                'CapSize', 18, 'Color', clr{i}, 'MarkerFaceColor', clr{i}, 'LineWidth', 1.5, ...
                'MarkerSize', markS(i), 'DisplayName', strrep(layer, '_', ': '));
    hold on
end

xticks(1:numel(params));
xticklabels(ticks);
xlabel('Parameter')
ylabel('Utilized (%)')
title({'Percentage of permanent structural traps utilized'})
legend()
drawnow

saveas(fig_util, strcat(plot_dir, 'struct_utilized'), 'png');

%% TESTING

%% Functions
function pc = runStandardPc(S, dummy_S, swr, snr, p_e, p_cap, layers, G)
    pc_vals = Capillary.PcNew(dummy_S, swr, snr, p_e, p_cap, 2);

    region_table = {[dummy_S, zeros(numel(dummy_S), 1)], [dummy_S, pc_vals]}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg(region_table, S, region_idx);
end

function pc = runLeverettJ_2D(S, dummy_S, phi, K, dummy_K, K_base, layers, G)
    [grid_Sw, grid_K] = ndgrid(dummy_S, dummy_K);
    pc_vals = Capillary.LeverettJ(grid_Sw, phi, grid_K, K_base);
    
    region_table = {{grid_Sw, grid_K, zeros(size(grid_Sw))}, ...
                     {grid_Sw, grid_K,  pc_vals}}; % container for pc values in each region
    region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
    pc = interpReg2D(region_table, S, K, region_idx);
end
