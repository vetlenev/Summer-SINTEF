%% Case 4 - gaussian relperm distribution & capillary pressure
% Simulates CO2 migration through multiple lowperm cells
% (more than case 1) on a large fine-scale grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Heterogeneous permeability - low k in cell layers, rest of domain high k.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui test-suite
ROOTDIR = strrep(ROOTDIR, '\', '/');
seed = rng();
seed.Seed = 5422; % Must set here, otherwise not updated

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
n_lowperm_layers = 30;
leaked_perc = 0.05; % allow 5% leakage

plot_base_dir = strcat(ROOTDIR, '../summer_sintef/case4/plots_noleakage');
data_dir = strcat(ROOTDIR, '../summer_sintef/case4/data_noleakage');
leaked_perc_str = erase(string(leaked_perc), '.');
plot_dir = sprintf(strcat(plot_base_dir, '/layers_%d_lp_%s'), n_lowperm_layers, leaked_perc_str);
dir_exists = mkdir(plot_base_dir) & mkdir(data_dir) & mkdir(plot_dir);

%% Define low-perm cells
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);

rand_x_start = linspace(perc_val(x, 0.0), perc_val(x, 0.9), n_lowperm_layers);
rand_x_start = rand_x_start(randperm(length(rand_x_start))); % random permute

rand_x_stop = rand_x_start + randi([round(perc_val(x, 0.1)), ...
                                    round(perc_val(x, 0.4))], ...
                                    [n_lowperm_layers, 1]).';
rand_x_stop = min(rand_x_stop, round(perc_val(x, 1))); % prevent layer going out of bounds

anticline_idx = 1:2:n_lowperm_layers; % every third layer is anticline

line_idx = setdiff(1:n_lowperm_layers, anticline_idx); % indices for line layers

z_start = linspace(perc_val(z, 0.05), perc_val(z, 0.9), n_lowperm_layers);
z_stop = z_start + perc_val(z, 0.02);

%% Generate layers
lowperm_cells = {};
all_lowperm_cells = [];
trapped_cells = {}; % structural trapping

corr_len_x = mean(rand_x_stop - rand_x_start) / 100;
corr_len_z = mean(z_stop - z_start) / 10;

% Make anticline forms
for i=anticline_idx
    num_x_cells = numel(G.cells.indexMap(x > rand_x_start(i) ... 
                        & x < rand_x_stop(i) ...
                        & z == min(z))); % only select ONE horizontal patch
    theta0 = pi/4;
    theta = linspace(theta0, pi - theta0, num_x_cells);
    r = (rand_x_stop(i) - rand_x_start(i))/2;
    z_anticline_start = -r*sin(theta) + z_start(i) + r/2*(1+cos(pi/2-theta0)); % -r*sin to get anticline (since z positive downwards)
    z_anticline_stop = -r*sin(theta) + z_stop(i) + r/2*(1+cos(pi/2-theta0)); % + r*cos(pi/4) to shift to bottom of spherical cap
    
    anticline_cells = [];
    anticline_cells_dummy = G.cells.indexMap(x > rand_x_start(i) & ... % minimum rectangle embedding entire curves layer
                                    x < rand_x_stop(i) & ...
                                    z > min(z_anticline_start) & ...
                                    z < max(z_anticline_stop));
    max_z = max(kk(anticline_cells_dummy));
    min_z = min(kk(anticline_cells_dummy));       
    z_bottom = unique(z(kk == max_z + 1)); % z-coord bottom part of layer
    
    trapped_cells{i} = [];
    
    for j=1:num_x_cells
        x_slice_start = (rand_x_stop(i) - rand_x_start(i))*(j-1)/num_x_cells + rand_x_start(i);
        x_slice_stop = (rand_x_stop(i) - rand_x_start(i))*j/num_x_cells + rand_x_start(i);
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
        trapped_cells{i} = cat(2, trapped_cells{i}, trapped_slice);
        
        anticline_cells = cat(1, anticline_cells, anticline_single);
    end

    num_z_cells = max_z-min_z+1; % +1 to get correct nr vertical cells in anticline form                                                  
    disp(numel(anticline_cells_dummy))
    disp(num_x_cells*num_z_cells)
    
    perm(anticline_cells_dummy) = lowperm + lowperm*FastGaussian([num_x_cells num_z_cells], 0.5, [corr_len_x corr_len_z]); % rectangular box around anticline
    perm(setdiff(anticline_cells_dummy, anticline_cells)) = baseperm; % reset perm for parts of rectangle NOT part of anticline
    all_lowperm_cells = cat(1, all_lowperm_cells, anticline_cells);
end

% Make straight forms
for i=line_idx
    if i == line_idx(fix(numel(line_idx)/2)) % put full-extended barrier in middle of domain
        rand_x_start(i) = min(x);
        rand_x_stop(i) = max(x);
        %z_stop(i) = z_start(i) + perc_val(z, 0.04); % increase vertical extent of the full-blocking layer
    end
    
    lowperm_cells{i} = G.cells.indexMap(x >= rand_x_start(i) & ...
                                    x <= rand_x_stop(i) & ...
                                    z > z_start(i) & ...
                                    z < z_stop(i));
    
    max_z = max(kk(lowperm_cells{i}));
    z_bottom = unique(z(kk == max_z + 1));
                                
    trapped_cells{i} = G.cells.indexMap(x >= rand_x_start(i) & ...
                                    x <= rand_x_stop(i) & ...
                                    z > z_stop(i) & ...
                                    z <= z_bottom).';   
                                
    [ix, jy, kz] = gridLogicalIndices(G, lowperm_cells{i});
    nxi = numel(unique(ix)); nzi = numel(unique(kz)); % dimensions of particular low-perm layer
    
    perm(lowperm_cells{i}) = lowperm + lowperm*FastGaussian([nxi nzi], 0.5, [corr_len_x corr_len_z]);
    all_lowperm_cells = cat(1, all_lowperm_cells, lowperm_cells{i});
end

% Make 

perm = max(perm, baseperm/1000); % perm must be positive -> cap at non-zero value to avoid singular matrix

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
p_e = 0.3*barsa;
p_cap = 3*barsa;
pc_vals = UtilFunctions.Pc(dummy_Sw, swr, p_e, p_cap, 2);

region_table = {[dummy_Sw, zeros(numel(dummy_Sw), 1)], [dummy_Sw, pc_vals]}; % container for pc values in each region
region_idx = {setdiff(G.cells.indexMap, all_lowperm_cells).', all_lowperm_cells}; % region to interpolate (rest, lowperm)
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
all_trapped_cells = cell2mat(trapped_cells);

clf;
f1 = UtilFunctions.fullsizeFig(1);
plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
perm_dummy = convertTo(perm, milli*darcy);
plotCellData(G, log10(perm_dummy), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'blue', 'EdgeColor', 'none');
colormap(autumn);
colorbarHist(log10(perm_dummy(all_lowperm_cells)), [min(log10(perm_dummy)), max(log10(perm_dummy))], 'South', 51);
title('Log of permeability field');
axis equal tight
view([0, 0])
zlim([min(z) max(z)]);
drawnow
hold off

saveas(f1, strcat(plot_dir, '/perm'), 'png');        

%% Plot structurally trapped cells
f2 = UtilFunctions.fullsizeFig(2);
plotGrid(G, all_lowperm_cells, 'FaceColor', 'yellow', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
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
top_cells = G.cells.indexMap(z < min(z(all_lowperm_cells)));
interior_cells = G.cells.indexMap(z >= min(z(all_lowperm_cells)));

p_top = fluid.rhoWS * norm(gravity) * min(z);
bc = pside(bc, G, 'Top', p_top, 'sat', [1 0]);

pz = fluid.rhoWS * norm(gravity) * unique(z); % hydrostatic pressure in entire domain
bc = pside(bc, G, 'Right', pz, 'sat', [1 0]);

tot_time = 8000*day();
dt = rampupTimesteps(tot_time, 70*day(), 10);

%% Initial state
% To simulate CO2 in supercritical phase, use initial pressure of 100 barsa
state = initResSol(G, 100*barsa, [1,0]);
t = 0;

f3 = UtilFunctions.fullsizeFig(3); % to hold saturations

plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, dummyW.cells, 'FaceColor', 'black', 'EdgeColor', 'none');
colormap(flipud(winter)); colorbar('southoutside'); caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, strcat(plot_dir, '/sat_0'), 'png');

f4 = UtilFunctions.fullsizeFig(4); % to hold cap pressure

plotGrid(G, all_lowperm_cells, 'FaceColor', 'none', 'EdgeColor', 'black', 'EdgeAlpha', 0.2);
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
pv = 100;
rates = [(sum(poreVolume(G, rock))/pv) / (inj_stop_rate*tot_time)];

[max_volumes, leaked_boundary, has_leaked, ...
    states, rates, rel_diff, categorized_vols] = FindMaxVolumeNoLeakage(G, rock, rates, state, model, ...
                                                                         all_trapped_cells, snr, dt, bc, ...
                                                                         inj_stop_rate, leaked_perc);


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
title({'Optimizing total injection volume.', sprintf('Green: below %.1f %% leakage.', leaked_perc*100), sprintf('Red: above %.1f %% leakage.', leaked_perc*100)});
drawnow
saveas(f5, strcat(plot_dir, '/opt_vol_search_', date), 'png');


%% Store volume categories
residual_ratio_filename = sprintf(strcat(data_dir, '/residual_leaked_%s_layers_%d_seed_%d.mat'), leaked_perc_str, n_lowperm_layers, seed.Seed);
struct_ratio_filename = sprintf(strcat(data_dir, '/struct_leaked_%s_layers_%d_seed_%d.mat'), leaked_perc_str, n_lowperm_layers, seed.Seed);
free_ratio_filename = sprintf(strcat(data_dir, '/free_leaked_%s_layers_%d_seed_%d.mat'), leaked_perc_str, n_lowperm_layers, seed.Seed);

residual_vol = categorized_vols{1};
structural_vol = categorized_vols{2};
free_vol = categorized_vols{3};

save(residual_ratio_filename, 'residual_vol'); % save to compare for different nr low-perm layers
save(struct_ratio_filename, 'structural_vol');
save(free_ratio_filename, 'free_vol');

