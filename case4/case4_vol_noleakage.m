%% Case 4 - gaussian relperm distribution & capillary pressure
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');
seed = rng();
seed.Seed = 6551; % Must set here, otherwise not updated

%% Define 2D grid
nx = 75; ny = 1; nz = 60; % 100 50
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

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/plots_optimal');
data_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/data_optimal');
leaked_perc_str = erase(string(leaked_perc), '.');
plot_dir = sprintf(strcat(plot_base_dir, '/lp_%s/layers_%d'), leaked_perc_str, n_lowperm_layers);
data_dir = sprintf(strcat(data_base_dir, '/lp_%s/layers_%d'), leaked_perc_str, n_lowperm_layers);
dir_exists = mkdir(plot_base_dir) & mkdir(data_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Define low-perm cells
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);

x_start = struct; x_stop = struct;
z_start = struct; z_stop = struct;

rand_x_start = linspace(perc_val(x, 0.0), perc_val(x, 0.9), n_lowperm_layers+n_imperm_layers);
rand_x_start = rand_x_start(randperm(length(rand_x_start))); % random permute

x_start.lowperm = rand_x_start(1:n_lowperm_layers);
x_start.imperm = rand_x_start(n_lowperm_layers+1:n_lowperm_layers+n_imperm_layers);

x_stop.lowperm = x_start.lowperm + randi([round(perc_val(x, 0.1)), ...
                                    round(perc_val(x, 0.4))], ...
                                    [n_lowperm_layers, 1]).';
x_stop.imperm = x_start.imperm + randi([round(perc_val(x, 0.1)), ...
                                    round(perc_val(x, 0.4))], ...
                                    [n_imperm_layers, 1]).';
                                
x_stop.lowperm = min(x_stop.lowperm, round(perc_val(x, 1))); % prevent layer going out of bounds
x_stop.imperm = min(x_stop.imperm, round(perc_val(x, 1)));

anticline_idx = 1:3:(n_lowperm_layers+n_imperm_layers); % every third layer is anticline

line_idx = setdiff(1:(n_lowperm_layers+n_imperm_layers), anticline_idx); % indices for line layers

rand_z_start = linspace(perc_val(z, 0.05), perc_val(z, 0.9), n_lowperm_layers+n_imperm_layers);
rand_z_start = rand_z_start(randperm(length(rand_z_start)));
z_start.lowperm = rand_z_start(1:n_lowperm_layers);
z_start.imperm = rand_z_start(n_lowperm_layers+1:n_lowperm_layers+n_imperm_layers);

z_stop.lowperm = z_start.lowperm + perc_val(z, 0.02);
z_stop.imperm = z_start.imperm + perc_val(z, 0.02);

%% Generate layers
added_layers = {[], []};
layer_types = {'lowperm', 'imperm'};
anticline_idxs = {anticline_idx(1:round(n_lowperm_layers/(n_lowperm_layers+n_imperm_layers)*length(anticline_idx))), ...
                    anticline_idx(round(n_lowperm_layers/(n_lowperm_layers+n_imperm_layers)*length(anticline_idx))+1:end) ...
                     - n_lowperm_layers};
line_idxs = {line_idx(1:round(n_lowperm_layers/(n_lowperm_layers+n_imperm_layers)*length(line_idx))), ...
                    line_idx(round(n_lowperm_layers/(n_lowperm_layers+n_imperm_layers)*length(line_idx))+1:end) ...
                     - n_lowperm_layers};                
perm_vals = {lowperm, 1e-5*milli*darcy};                

trapped_cells = {}; % structural trapping
corr_len_x = mean(x_stop.lowperm - x_start.lowperm) / 100;
corr_len_z = mean(z_stop.lowperm - z_start.lowperm) / 10;

% Make anticline forms
for k=1:numel(layer_types)
    layer_type = layer_types{k};
    for i=anticline_idxs{k}
        num_x_cells = numel(G.cells.indexMap(x > x_start.(layer_type)(i) ... 
                            & x < x_stop.(layer_type)(i) ...
                            & z == min(z))); % only select ONE horizontal patch
        theta0 = pi/4;
        theta = linspace(theta0, pi - theta0, num_x_cells);
        r = (x_stop.(layer_type)(i) - x_start.(layer_type)(i))/2;
        z_anticline_start = -r*sin(theta) + z_start.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % -r*sin to get anticline (since z positive downwards)
        z_anticline_stop = -r*sin(theta) + z_stop.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % + r*cos(pi/4) to shift to bottom of spherical cap

        anticline_cells = [];
        anticline_cells_dummy = G.cells.indexMap(x > x_start.(layer_type)(i) & ...
                                        x < x_stop.(layer_type)(i) & ...
                                        z > min(z_anticline_start) & ...
                                        z < max(z_anticline_stop));
                                    
        max_z = max(kk(anticline_cells_dummy));
        min_z = min(kk(anticline_cells_dummy)); 
        z_bottom = unique(z(kk == max_z + 1)); % z-coord bottom part of layer
        
        trapped_cells_cat = [];
        
        for j=1:num_x_cells
            x_slice_start = (x_stop.(layer_type)(i) - x_start.(layer_type)(i))*(j-1)/num_x_cells + x_start.(layer_type)(i);
            x_slice_stop = (x_stop.(layer_type)(i) - x_start.(layer_type)(i))*j/num_x_cells + x_start.(layer_type)(i);
            anticline_single = G.cells.indexMap(x > x_slice_start & ...
                                        x < x_slice_stop & ...
                                        z > z_anticline_start(j) & ...
                                        z < z_anticline_stop(j));
                                    
            [ix, jy, kz] = gridLogicalIndices(G, anticline_single);       
            %max_z = max(max(kz), max_z); % max depth of this anticline form
            %min_z = min(min(kz), min_z); % min depth of this anticline form
            trapped_slice = G.cells.indexMap(x > x_slice_start & ...
                                         x < x_slice_stop & ...
                                         z > z_anticline_stop(j) & ...
                                         z <= z_bottom).';
            trapped_cells_cat = cat(2, trapped_cells_cat, trapped_slice);

            anticline_cells = cat(1, anticline_cells, anticline_single);
        end
        
        trapped_cells = cat(1, trapped_cells, trapped_cells_cat);
        
        num_z_cells = max_z-min_z+1; % +1 to get correct nr vertical cells in anticline form            

        overlapped_cells = ismember(anticline_cells_dummy, vertcat(added_layers{:}));
        store_perm = perm(overlapped_cells);
        
        %anticline_cells_unique = setdiff(anticline_cells_dummy, vertcat(added_layers{:}));
        perm(anticline_cells_dummy) = perm_vals{k} + perm_vals{k}*FastGaussian([num_x_cells num_z_cells], 0.8, [corr_len_x corr_len_z]); % rectangular box around anticline
        perm(overlapped_cells) = store_perm;
        perm(setdiff(anticline_cells_dummy, vertcat(anticline_cells, vertcat(added_layers{:})))) = baseperm; % reset perm for parts of rectangle NOT part of anticline
        added_layers{k} = cat(1, added_layers{k}, anticline_cells);
    end
    
    % Make straight forms
    for i=line_idxs{k}
        if i == line_idxs{k}(fix(end/2)) && strcmp(layer_type, 'lowperm')
            x_start.(layer_type)(i) = min(x);
            x_stop.(layer_type)(i) = max(x);
            z_stop.(layer_type)(i) = z_start.(layer_type)(i) + perc_val(z, 0.04); % increase vertical extent of the full-blocking layer
        end

        line_cells = G.cells.indexMap(x >= x_start.(layer_type)(i) & ...
                                        x <= x_stop.(layer_type)(i) & ...
                                        z > z_start.(layer_type)(i) & ...
                                        z < z_stop.(layer_type)(i));
                                    
        max_z = max(kk(line_cells));
        z_bottom = unique(z(kk == max_z + 1));
        
        trapped_cells_cat = G.cells.indexMap(x >= x_start.(layer_type)(i) & ...
                                    x <= x_stop.(layer_type)(i) & ...
                                    z > z_stop.(layer_type)(i) & ...
                                    z <= z_bottom).';
                                
        trapped_cells = cat(1, trapped_cells, trapped_cells_cat);
                                    
        [ix, jy, kz] = gridLogicalIndices(G, line_cells);
        nxi = numel(unique(ix)); nzi = numel(unique(kz)); % dimensions of particular low-perm layer
        
        perm(line_cells) = perm_vals{k} + perm_vals{k}*FastGaussian([nxi nzi], 0.8, [corr_len_x corr_len_z]);
        added_layers{k} = cat(1, added_layers{k}, line_cells);
    end
end

all_added_layers = vertcat(added_layers{:});
perm(all_added_layers) = max(perm(all_added_layers), baseperm/100); % perm must be positive -> cap at non-zero value to avoid singular matrix

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
pc_vals = UtilFunctions.Pc(dummy_Sw, swr, p_e, p_cap, 2);

region_table = {[dummy_Sw, zeros(numel(dummy_Sw), 1)], [dummy_Sw, pc_vals]}; % container for pc values in each region
region_idx = {setdiff(G.cells.indexMap, all_added_layers).', all_added_layers}; % region to interpolate (rest, lowperm)
fluid.pcOW = @(S, varargin) interpReg(region_table, S, region_idx); % in both regions, interpolate over saturations S
% or using Leverett-J pc to account for differences in permeability

figure(1);
plot(dummy_Sw, pc_vals, 'LineWidth', 1.5);
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

%% Boundary conditions, schedule
bc = []; % no-flux as default
top_cells = G.cells.indexMap(z < min(z(all_added_layers)));
interior_cells = G.cells.indexMap(z >= min(z(all_added_layers)));

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, G, 'Right', pz, 'sat', [1 0]);

tot_time = 16000*day();
dt = rampupTimesteps(tot_time, 75*day(), 10);

inj_years = regexp(formatTimeRange(tot_time), '\d+ Years', 'match');
years = strrep(inj_years, ' ', '_');

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
colorbar; caxis([0, 2*p_e]);
title({'Capillary pressure (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f4, strcat(plot_dir, '/cap_pres_0'), 'png');

%% Run experiment
inj_stop_rate = 0.3;
% start with injection rate corresponding to 1/pv'th of total pore volume
pv = 100; % start with rate yielding total volume of 1/100th of pore volume in domain
rates = [(sum(poreVolume(G, rock))/pv) / (inj_stop_rate*tot_time)];
kmax = 12; % max 10 optimization steps

[max_volumes, leaked_boundary, has_leaked, ...
    states, rates, scales, rel_diff, categorized_vols] = FindMaxVolumeNoLeakage(G, rock, rates, state, model, ...
                                                                         all_trapped_cells, snr, dt, bc, ...
                                                                         inj_stop_rate, leaked_perc, kmax);


%% Show optimal solution found                                                       
disp(max_volumes)
dt_plot = cat(1, repmat(5, [fix(numel(states)/4), 1]), ...
                 repmat(10, [numel(states)-fix(numel(states)/4), 1]));

for i=1:numel(states)
    t = t + dt(i);   

    if ~mod(i, 10)
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
        colorbar; caxis([0, 2*p_e]);
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
title({'Optimizing total CO2 volume, injecting for', formatTimeRange(tot_time), sprintf('Green: below %.1f %% leakage. Red: above %.1f %% leakage.', leaked_perc*100, leaked_perc*100)});
drawnow
saveas(f5, strcat(plot_dir, '/opt_vol_search_', date), 'png');


%% Store volume categories
%residual_ratio_filename = sprintf(strcat(data_dir, '/residual_leaked_%s_layers_%d_seed_%d.mat'), leaked_perc_str, n_lowperm_layers, seed.Seed);
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_%s_seed_%d.mat'), years{1}, seed.Seed);
struct_ratio_filename = sprintf(strcat(data_dir, '/struct_%s_seed_%d.mat'), years{1}, seed.Seed);
free_ratio_filename = sprintf(strcat(data_dir, '/free_%s_seed_%d.mat'), years{1}, seed.Seed);
inj_rate_filename = sprintf(strcat(data_dir, '/inj_rate_%s_seed_%d.mat'), years{1}, seed.Seed);

residual_vol = categorized_vols{1};
structural_vol = categorized_vols{2};
free_vol = categorized_vols{3};
inj_rate = rates(end);

save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
save(struct_ratio_filename, 'structural_vol');
save(free_ratio_filename, 'free_vol');
save(inj_rate_filename, 'inj_rate');

%% Read data and store in structs
lp_folders = dir(strcat(data_base_dir, '/lp_*'));
lp_folders = UtilFunctions.sortStructByField(lp_folders, 'name');

num_lp = numel(lp_folders);
num_layers = zeros(num_lp, 1);

used_lab = {{}, {{}}}; % first: percentage leakage, second: num layers
unique_lab_idx = {[], []};

structural_struct = struct;
residual_struct = struct;
free_struct = struct;
inj_struct = struct;

for i=1:num_lp
    lab_leaked = regexp(lp_folders(i).name, 'lp_\d+', 'match');     
    used_lab{1} = cat(2, used_lab{1}, lab_leaked{1});
    %unique_lab_idx{1} = cat(2, unique_lab_idx{1}, i);
    
    %if lab_leaked{1} == leaked_perc_str % only plot layer configuration for current leaked percentage               
    layer_folders = dir(strcat(lp_folders(i).folder, '/', lp_folders(i).name, '/layers_*'));  
    layer_folders = UtilFunctions.sortStructByField(layer_folders, 'name');
    num_layers(i) = numel(layer_folders);   

    for j=1:num_layers(i)
        lab_layers = regexp(layer_folders(j).name, 'layers_\d+', 'match');
        %if ~any(ismember(used_lab{2}{i}, lab_layers{1}))
        used_lab{2}{i} = cat(2, used_lab{2}{i}, lab_layers{1});      
        
        inj_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/inj_rate_*.mat'));
        inj_struct.(lab_leaked{1}).(lab_layers{1}) = [];
        
        for k=1:numel(inj_files)
            load_injrate = load(strcat(inj_files(k).folder, '\', inj_files(k).name), 'inj_rate');
            inj_struct.(lab_leaked{1}).(lab_layers{1}) = cat(1, inj_struct.(lab_leaked{1}).(lab_layers{1}), load_injrate.inj_rate);
        end

        structural_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/struct_*.mat'));
        free_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/free_*.mat'));
        residual_files = dir(strcat(layer_folders(j).folder, '/', layer_folders(j).name, '/residual_*.mat'));
        
        structural_struct.(lab_leaked{1}).(lab_layers{1}) = [];
        free_struct.(lab_leaked{1}).(lab_layers{1}) = [];
        residual_struct.(lab_leaked{1}).(lab_layers{1}) = [];
        
        for k=1:numel(free_files) % assume same number of files for given lp and num layers
            load_structural = load(strcat(structural_files(k).folder, '\', structural_files(k).name), 'structural_vol');
            % NB: appended row-wise, so to access a specific seed do:
            % structural_struct.layers_X.lp_Y(i,:)
            structural_struct.(lab_leaked{1}).(lab_layers{1}) = cat(2, structural_struct.(lab_leaked{1}).(lab_layers{1}), load_structural.structural_vol);       
            load_free = load(strcat(free_files(k).folder, '\', free_files(k).name), 'free_vol');
            free_struct.(lab_leaked{1}).(lab_layers{1}) = cat(2, free_struct.(lab_leaked{1}).(lab_layers{1}), load_free.free_vol);  
            load_residual = load(strcat(residual_files(k).folder, '\', residual_files(k).name), 'residual_vol');
            residual_struct.(lab_leaked{1}).(lab_layers{1}) = cat(2, residual_struct.(lab_leaked{1}).(lab_layers{1}), load_residual.residual_vol);  
        end
    end
end


%% Plot optimal injection rate across lp's
% Compute mean and variance for each num layers
mean_lp = struct;
var_lp = struct;
ticks = {};

for i=1:num_lp
    lp_lab = used_lab{1}(i);
    ticks{i} = strrep(lp_lab{1}, 'lp_0', '0.');
    for j=1:num_layers(i)
        layer_lab = used_lab{2}{i}(j);
        mean_lp.(layer_lab{1}).(lp_lab{1}) = mean(inj_struct.(lp_lab{1}).(layer_lab{1}), 1);
        var_lp.(layer_lab{1}).(lp_lab{1}) = var(inj_struct.(lp_lab{1}).(layer_lab{1}), 0, 1);

    end
end


f6 = figure(6);
fields_layers = fieldnames(mean_lp);
for i=1:numel(fields_layers)
    plot(1:numel(mean_lp.(fields_layers{i})), struct2array(mean_lp.(fields_layers{i})), ...
        '--.', 'MarkerSize', 15, 'DisplayName', strrep(fields_layers{i}, '_', ' '));
    hold on    
end

xlabel('Allowed leakage (percentage)')
ylabel('Rate (m^3/s)')
xticks(1:num_lp);
xticklabels(ticks);

tot_time_cut = regexp(formatTimeRange(tot_time), '.+?Days', 'match');
inj_time_cut = regexp(formatTimeRange(tot_time*inj_stop_rate), '.+?Days', 'match');
title({'Optimal injection rate to satifsy given leakage amount.', ...
        strcat('Simulation time: ', tot_time_cut{1}), ...
        strcat('Injection time: ', inj_time_cut{1})})
legend();    
drawnow

saveas(f6, strcat(plot_base_dir, '/all_optimal_rates'), 'png');
    

%% Plot volume distributions
fm = 7;
fig_mean = figure(fm);
fv = 8;
fig_var = figure(fv);

linS = {'-', '--', ':', '-o', '--o'};

for i=1:num_lp % one figure for each lp
    mean_layers = struct;
    var_layers = struct;
    
    lp_lab = used_lab{1}(i);
    lp_perc = strrep(lp_lab{1}, 'lp_0', '0.');
    for j=1:num_layers(i)
        layer_lab = used_lab{2}{i}(j);
        mean_layers.structural.(layer_lab{1}) = mean(structural_struct.(lp_lab{1}).(layer_lab{1}), 2);
        var_layers.structural.(layer_lab{1}) = var(structural_struct.(lp_lab{1}).(layer_lab{1}), 0, 2);
        mean_layers.free.(layer_lab{1}) = mean(free_struct.(lp_lab{1}).(layer_lab{1}), 2);
        var_layers.free.(layer_lab{1}) = var(free_struct.(lp_lab{1}).(layer_lab{1}), 0, 2);
        mean_layers.residual.(layer_lab{1}) = mean(residual_struct.(lp_lab{1}).(layer_lab{1}), 2);
        var_layers.residual.(layer_lab{1}) = var(residual_struct.(lp_lab{1}).(layer_lab{1}), 0, 2);
        
        %t = 1:numel(mean_layers.structural.(layer_lab{1}));
        figure(fm);
        hold on      
        plot(mean_layers.structural.(layer_lab{1}), linS{j}, 'Color', 'blue', 'DisplayName', ['Structural: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5);      
        plot(mean_layers.residual.(layer_lab{1}), linS{j}, 'Color', 'red', 'DisplayName', ['Residual: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5); 
        plot(mean_layers.free.(layer_lab{1}), linS{j}, 'Color', 'green', 'DisplayName', ['Free: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5); 
        
        figure(fv);
        hold on
        plot(var_layers.structural.(layer_lab{1}), 'Color', 'blue', 'DisplayName', ['Structural: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5);
        plot(var_layers.residual.(layer_lab{1}), 'Color', 'red', 'DisplayName', ['Residual: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5);
        plot(var_layers.free.(layer_lab{1}), 'Color', 'green', 'DisplayName', ['Free: ', strrep(layer_lab{1}, '_', ' ')], 'LineWidth', 1.5);
               
    end
    
    figure(fm)
    xlabel('Time step')
    ylabel('Volume (m^3)')
    title({'Mean volume distribution.', ['Allowing ', lp_perc, ' rate of leakage.']})  
    legend('Location', 'northwest')  
    saveas(fig_mean, sprintf(strcat(plot_base_dir, '/lp_%s/mean_trapping'), leaked_perc_str), 'png');
    
    figure(fv)
    xlabel('Time step')
    ylabel('Volume (m^3)')
    title({'Variance of volume distribution.', ['Allowing ', lp_perc, ' rate of leakage.']}) 
    legend()
    saveas(fig_var, sprintf(strcat(plot_base_dir, '/lp_%s/var_trapping'), leaked_perc_str), 'png')
    
    fm = fm + 2;
    fig_mean = figure(fm);
    fv = fv + 2;
    fig_var = figure(fv);
end


%% TESTING
