%% Case 4 - gaussian relperm distribution & capillary pressure
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');

my_seed = 6316;
assert_timestep = true;

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
lowperm = 20*milli*darcy;
baseperm = 100*milli*darcy;
perm = repmat(baseperm, [G.cells.num 1]);
poro = 0.3;

%% Directories
n_lowperm_layers = 15;
n_imperm_layers = round(n_lowperm_layers/3);
leaked_perc = 0.0; % allow X*100% leakage

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/plots_timestep');
data_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/data_timestep');
leaked_perc_str = erase(string(leaked_perc), '.');
plot_dir = sprintf(strcat(plot_base_dir, '/lp_%s/layers_%d'), leaked_perc_str, n_lowperm_layers+n_imperm_layers);
data_dir = sprintf(strcat(data_base_dir, '/lp_%s/layers_%d'), leaked_perc_str, n_lowperm_layers+n_imperm_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

this_lp = regexp(data_dir, 'lp_\d+', 'match');
this_layer = regexp(data_dir, 'layers_\d+', 'match');

%% Set seed
seed = UtilFunctions.setSeed(data_dir, my_seed);
rng(seed)

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
hold off
                       
%% Capillary pressure
dummy_Sw = linspace(0, 1, G.cells.num)';
p_e = 0.5*barsa;
p_cap = 3*barsa;

median_pc = 2*barsa;
%S_scaled = max((1-dummy_Sw-swr)./(1-swr), 1e-5);
standard_pc = 0;

if standard_pc
    plot_pc_lim = p_cap;
    fluid.pcOW = @(S, varargin) Capillary.runStandardPc(S, dummy_Sw, swr, snr, p_e, p_cap, all_added_layers, G);
else
    plot_pc_lim = 5*median_pc;
    fluid.pcOW = @(S, varargin) Capillary.LeverettJ(S, poro, perm, baseperm, median_pc);
end

figure(1);
plot(dummy_Sw, fluid.pcOW(dummy_Sw), 'LineWidth', 1.5);
xlabel('Water saturation');
title('Capillary pressure function');
saveas(f1, strcat(plot_dir, '/cap_pres'), 'png');
hold off

%% Dummy well
well_h = 1; % cell perforations in vertical direction
perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
dummyW = addWell([], G, rock, perforation_idx, ...
        'Type', 'rate', 'Val', 1*meter^3/day(), ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');
    
%% Plot permeability field
all_trapped_cells = horzcat(trapped_cells{:});

clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(perm, milli*darcy);
plotCellData(G, log10(perm_dummy), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn);
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
plotGrid(G, all_added_layers, 'FaceColor', 'yellow', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotGrid(G, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
plotGrid(G, all_trapped_cells, 'FaceColor', [0.5, 0, 0.5], 'EdgeColor', 'none');
title('Structurally trapped CO2');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
drawnow
hold off

saveas(f2, strcat(plot_dir, '/struct_trapped'), 'png');  


%% Set up solver
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
disp(model)

%% Boundary conditions, schedule and initial conditions
top_cells = G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = G.cells.indexMap(z >= min(z(all_added_layers)));

bc = []; % no-flux as default

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, G, 'Right', pz, 'sat', [1 0]);

tot_time = 1000*365*day();
%dt = rampupTimesteps(tot_time, 1000*day(), 10);
dt = rampupTimesteps(round(tot_time/4), 300*day(), 10);
dt = [dt', rampupTimesteps(tot_time-round(tot_time/4), 800*day(), 0)']';
num_dt = numel(dt);

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

inj_stop_rate = 0.05;
rates = 0.64*meter^3/day();
kmax = 1; % only do one iteration - we are interested in comparison, not optimal

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(G, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

f4 = UtilFunctions.fullsizeFig(4); % to hold cap pressure

plotGrid(G, all_added_layers, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, fluid.pcOW(state.s(:,1)), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colorbar; caxis([0, plot_pc_lim]);
title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f4, strcat(plot_dir, '/cap_pres_0'), 'png');

%% Run experiment

[max_volumes, leaked_boundary, has_leaked, ...
    states, rates, scales, rel_diff, categorized_vols] = FindMaxVolumeNoLeakage(G, rock, rates, state, model, ...
                                                                         all_trapped_cells, swr, snr, dt, bc, ...
                                                                         inj_stop_rate, leaked_perc, kmax);


%% Show optimal solution found                                                       
disp(max_volumes)
dt_plot = cat(1, repmat(10, [fix(numel(states)/10), 1]), ...
                 repmat(30, [numel(states)-fix(numel(states)/10), 1]));

for i=1:numel(states)
    t = t + dt(i);   

    if ~mod(i, dt_plot(i))
        figure(3)
        set(f3, 'visible', 'off');
        plotCellData(G, states{i}.s(:,1), 'EdgeColor', 'none');
        plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
        colormap(flipud(winter)); colorbar('southoutside'); caxis([0 1]);
        title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
        axis equal tight;
        view([0, 0]);   

        set(f3, 'visible', 'on');
        filename_f3 = sprintf(strcat(plot_dir, '/sat_%d'), i);
        saveas(f3, filename_f3, 'png');
        
        figure(4);
        set(f4, 'visible', 'off');
        plotCellData(G, fluid.pcOW(states{i}.s(:,1)), 'EdgeColor', 'none');
        plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none'); 
        colorbar; caxis([0, plot_pc_lim]);
        title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
        axis equal tight
        view([0, 0])
       
        set(f4, 'visible', 'on');
        filename_f4 = sprintf(strcat(plot_dir, '/cap_pres_%d'), i);
        saveas(f4, filename_f4, 'png');
    end 
end

disp(leaked_boundary)

% Plot the search for max volume
f5 = figure(5);

%plot_rates = rates(1:numel(max_volumes))/(meter^3)*day();
plot_rates = 1:numel(max_volumes);
for ii=1:numel(max_volumes)
    ms = 25;
    if ii == numel(max_volumes) % final optimal volume
        c = 'p';
        ms = 12;
    elseif has_leaked(ii)
        c = '.r';
        lab = sprintf('above %d %% leakage', leaked_perc*100);
    else
        c = '.g';
    end
    plot(plot_rates(ii), max_volumes(ii), c, 'MarkerSize', ms, 'MarkerFaceColor', 'green');
    hold on
end

plot(plot_rates, max_volumes, 'black');

xlabel('Iterations');
ylabel('Volume (m^3)');
xlim([0, plot_rates(end)+1]);
title({'Optimizing total CO2 volume, injecting for', formatTimeRange(tot_time), sprintf('Green: below %.1f %% leakage. Red: above %.1f %% leakage.', leaked_perc*100, leaked_perc*100)});
drawnow
saveas(f5, strcat(plot_dir, '/opt_vol_search_', date), 'png');


%% Store volume categories
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_%s_seed_%d_dt_%d.mat'), years{1}, seed.Seed, num_dt);
struct_ratio_filename = sprintf(strcat(data_dir, '/struct_%s_seed_%d_dt_%d.mat'), years{1}, seed.Seed, num_dt);
free_ratio_filename = sprintf(strcat(data_dir, '/free_%s_seed_%d_dt_%d.mat'), years{1}, seed.Seed, num_dt);
inj_rate_filename = sprintf(strcat(data_dir, '/inj_rate_%s_seed_%d_dt_%d.mat'), years{1}, seed.Seed, num_dt);

residual_vol = categorized_vols{1};
structural_vol = categorized_vols{2};
free_vol = categorized_vols{3};
inj_rate = rates(end);

save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
save(struct_ratio_filename, 'structural_vol');
save(free_ratio_filename, 'free_vol');
save(inj_rate_filename, 'inj_rate');

dt_filename = sprintf(strcat(data_dir, '/dt_%d.mat'), num_dt);
save(dt_filename, 'dt');

%% Read data and store in structs
lp_folders = dir(strcat(data_base_dir, '/lp_*'));
lp_folders = UtilFunctions.sortStructByField(lp_folders, 'name');

num_lp = numel(lp_folders);
num_layers = zeros(num_lp, 1);

used_lab = {{}, {}, {}}; % first: percentage leakage, second: num layers, third: timestep
unique_lab_idx = {[], []};

structural_struct = struct; % struct to hold structurally trapped values for each leaked percentage and num layer
residual_struct = struct;
free_struct = struct;
inj_struct = struct;

for i=1:num_lp
    lab_leaked = regexp(lp_folders(i).name, 'lp_\d+', 'match');     
    used_lab{1} = cat(2, used_lab{1}, lab_leaked{1});
    used_lab{2}{i} = {};
    
    %if lab_leaked{1} == leaked_perc_str % only plot layer configuration for current leaked percentage               
    layer_folders = dir(strcat(lp_folders(i).folder, '/', lp_folders(i).name, '/layers_*'));  
    layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
    num_layers(i) = numel(layer_folders);   

    for j=1:num_layers(i)
        lab_layers = regexp(layer_folders(j).name, 'layers_\d+', 'match'); 
        used_lab{2}{i} = cat(1, used_lab{2}{i}, lab_layers{1}); 
        used_lab{3}{i}{j} = {};
        
        inj_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/inj_rate_*.mat'));
        inj_struct.(lab_leaked{1}).(lab_layers{1}) = [];
        
        structural_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/struct_*.mat'));
        free_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/free_*.mat'));
        residual_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/residual_*.mat'));
        
        structural_struct.(lab_leaked{1}).(lab_layers{1}) = {};
        free_struct.(lab_leaked{1}).(lab_layers{1}) = {};
        residual_struct.(lab_leaked{1}).(lab_layers{1}) = {};       
              
        for k=1:numel(inj_files)
            load_injrate = load(strcat(inj_files(k).folder, '\', inj_files(k).name), 'inj_rate');
            inj_struct.(lab_leaked{1}).(lab_layers{1}) = cat(1, inj_struct.(lab_leaked{1}).(lab_layers{1}), load_injrate.inj_rate * days/seconds); % NB: units of m^3/day
            
            free_file = strcat(free_files(k).folder, '\', free_files(k).name);
            lab_dt = regexp(free_file, 'dt_\d+', 'match');
            used_lab{3}{i}{j} = cat(1, used_lab{3}{i}{j}, lab_dt{1}); 
            
            load_structural = load(strcat(structural_files(k).folder, '\', structural_files(k).name), 'structural_vol');            
            structural_struct.(lab_leaked{1}).(lab_layers{1}).(lab_dt{1}) = load_structural.structural_vol;       
            load_free = load(strcat(free_files(k).folder, '\', free_files(k).name), 'free_vol');
            free_struct.(lab_leaked{1}).(lab_layers{1}).(lab_dt{1}) = load_free.free_vol;  
            load_residual = load(strcat(residual_files(k).folder, '\', residual_files(k).name), 'residual_vol');
            residual_struct.(lab_leaked{1}).(lab_layers{1}).(lab_dt{1}) = load_residual.residual_vol;
            % NB: appended row-wise, so to access a specific seed do:
            % structural_struct.layers_X.lp_Y(i,:)
        end
               
    end
end


%% Plot optimal injection rate across lp's
% Compute mean and variance for each num layers and leaked percentage
mean_lp = struct;
std_lp = struct;
ticks = {};

for i=1:num_lp
    lp_lab = used_lab{1}(i);   
    ticks{i} = strrep(lp_lab{1}, 'lp_0', '0.');  
    for j=1:num_layers(i)
        layer_lab = used_lab{2}{i}(j);
        mean_lp.(layer_lab{1}).(lp_lab{1}) = mean(inj_struct.(lp_lab{1}).(layer_lab{1}), 1);
        std_lp.(layer_lab{1}).(lp_lab{1}) = std(inj_struct.(lp_lab{1}).(layer_lab{1}), 0, 1);
        % to access index of current lp and num_layer - needed for plotting
        if strcmp(lp_lab{1}, this_lp{1}) && strcmp(layer_lab{1}, this_layer{1})
           this_lp_idx = i;
           this_layer_idx = j; 
        end
    end
end


f6 = figure(6);
fields_layers = fieldnames(mean_lp); % number of unique layers (across all leaked percentages)
clr = {'blue', 'red', 'green', 'magenta', 'orange'};

for i=1:numel(fields_layers)
    mean_lp_layer = mean_lp.(fields_layers{i});
    std_lp_layer = std_lp.(fields_layers{i});
    num_lp_layer = numel(fieldnames(mean_lp_layer)); % number of unique leakage percentages for this layer

    errorbar(1:num_lp_layer, struct2array(mean_lp_layer), struct2array(std_lp_layer), '-o', ...
            'Color', clr{i}, 'MarkerFaceColor', clr{i}, ...
            'DisplayName', ['layers: ', strrep(fields_layers{i}, 'layers_', '')]);
    ylabel('Mean (m^3/day)');
    hold on
end

xlabel('Allowed leakage (ratio)')
xticks(1:num_lp);
xticklabels(ticks);

tot_time_cut = regexp(formatTimeRange(tot_time), '.+?Days', 'match');
inj_time_cut = regexp(formatTimeRange(tot_time*inj_stop_rate), '.+?Days', 'match');
title({'Optimal injection rate to satifsy given leakage amount.', ...
        ['Simulation time: ', tot_time_cut{1}], ...
        ['Injection time: ', inj_time_cut{1}]})
legend('Location', 'northwest');
drawnow

saveas(f6, strcat(plot_base_dir, '/all_optimal_rates'), 'png');


%% Plot variation for different time step selection
% test to what degree results for long vs. short time steps coincide   
    
this_structural = structural_struct.(this_lp{1}).(this_layer{1});
this_residual = residual_struct.(this_lp{1}).(this_layer{1});
this_free = free_struct.(this_lp{1}).(this_layer{1});

fig_time = figure(7);
lp_perc = strrep(this_lp{1}, 'lp_0', '0.');
linS = {'-', '--', ':', '-o', '--o'};

for i=1:numel(fieldnames(this_free))
    this_dt = used_lab{3}{this_lp_idx}{this_layer_idx}(i);       
    this_dt_int = str2num(cell2mat(strrep(this_dt, 'dt_', '')));   
    this_dts_filename = sprintf(strcat(data_dir, '/dt_%d.mat'), this_dt_int);      
    this_dt_str = strrep(this_dt{1}, 'dt_', 'dt ');
    
    this_dts = load(this_dts_filename, 'dt');
    this_dts = cumsum(this_dts.dt)*second()/year();    
    this_dts = [0, this_dts']; 

    plot(this_dts, this_structural.(this_dt{1}), linS{i}, 'Color', 'blue', 'DisplayName', ['S: ', this_dt_str], 'LineWidth', 1.5);
    hold on
    plot(this_dts, this_residual.(this_dt{1}), linS{i}, 'Color', 'red', 'DisplayName', ['R: ', this_dt_str], 'LineWidth', 1.5);
    hold on
    plot(this_dts, this_free.(this_dt{1}), linS{i}, 'Color', 'green', 'DisplayName', ['F: ', this_dt_str], 'LineWidth', 1.5);
    hold on
end

legend();
hold off
xlabel('Time (years)');
ylabel('Volume (m^3)');
title({'Comparison of different number of time steps.', ...
        ['Layers: ', sprintf('%d', n_lowperm_layers+n_imperm_layers)], ...
        ['Allowed leakage rate: ', lp_perc]});   

saveas(fig_time, strcat(plot_dir, '/assert_timestep'), 'png');

%% TESTING
