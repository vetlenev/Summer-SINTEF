%% 2D buoyancy-driven CO2 migration - impermeable cells
% Simulates CO2 migration through multiple impermeable cells
% embedded in a rectangular grid.
% Rightmost boundary is open for water and CO2 to exit, others closed.
% CO2 injected at constant rate in lower left corner.
% Homogeneous permeability.
mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui

%% Define 2D grid
nx = 50; ny = 1; nz = 35;
lx = 80; ly = 1; lz = 60;
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

%% Add impermeable cells
seed = rng(2022);
n_imperm_cells = 6;
perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);

rand_x_start = randi([round(perc_val(x, 0.1)), round(perc_val(x, 0.7))], [n_imperm_cells,1]);
rand_x_stop = rand_x_start + randi([round(0.1*(max(x)-min(x))), ...
                                    round(0.4*(max(x)-min(x)))], ...
                                    [n_imperm_cells, 1]); %randi([rand_x_start + perc_val(x,0.1), perc_val(x,1)]);
rand_x_stop = min(rand_x_stop, round(perc_val(x, 1)));

rand_z = z(randi([0.1*length(z), 0.9*length(z)], [n_imperm_cells,1]));

imperm_cells = {};
all_imperm_cells = [];
for i=1:n_imperm_cells
    imperm_cells{i} = G.cells.indexMap(x > rand_x_start(i) & ...
                                    x < rand_x_stop(i) & ...
                                    z == rand_z(i));
    perm(imperm_cells{i}) = 1e-3*milli*darcy;
    all_imperm_cells = cat(1, all_imperm_cells, imperm_cells{i});
end

%% Compute rock+fluid objects
rock = makeRock(G, perm, poro);
T = computeTrans(G, rock, 'Verbose', true);

fluid = initSimpleADIFluid('phases', 'WO', ... % [water, GAS] or OIL?
                           'mu', [1, 0.015]*centi*poise, ... % viscosity
                           'n',  [2, 2], ... % relperm powers
                           'rho', [1000, 2]*kilogram/meter^3); % densities: [water, CO2]                                                          

%% Horizontal well
rate = 7*meter^3/day; % 7
% Put well slightly above bottom to avoid it interacting with bottom BC
perforation_idx = G.cells.indexMap(z == max(z)-lz/nz & x < 10);
W = addWell([], G, rock, perforation_idx, ...
            'Type', 'rate', 'Val', rate, ...
            'InnerProduct', 'ip_tpf', ...
            'Radius', 0.1, 'Dir', 'x', ...
            'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
            'Name', 'P1');

%W.lims.bhp = 1*barsa;
%% Plot grid
% Map to associated faces
imperm_faces = zeros([size(all_imperm_cells), 6]); % 6 faces in 3d cartesian coords
for i=1:size(all_imperm_cells)
    icell_start = G.cells.facePos(all_imperm_cells(i));
    icell_stop = G.cells.facePos(all_imperm_cells(i)+1)-1;
    imperm_faces(i,:) = G.cells.faces(icell_start:icell_stop, 1);
end

f1 = figure(1);
plotFaces(G, imperm_faces(:), 'FaceColor', 'none', 'EdgeColor', 'black', 'Linewidth', 1.5);
plotCellData(G, perm);
plotGrid(G, W.cells, 'FaceColor', 'blue', 'EdgeColor', 'k', 'LineWidth', 1);
colormap(autumn); colorbar;
title('Permeability field (milli*darcy)');
axis equal tight
view([0, 0])
drawnow
hold off

saveas(f1, 'summer_sintef/case1/plots/perm_field', 'png');

%% Set up solver
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
disp(model)

%% Boundary conditions and schedule
bc = []; % no-flux as default
%bc = fluxside(bc, G, 'EAST', -10*meter^3/day, 'sat', [0 1]); % open right boundary (minus for fluc to be OUT)
pz = fluid.rhoWS * norm(gravity) * uniquetol(z); % hydrostatic pressure (uniquetol necessary to avoid some duplicates by round-off)
bc = pside(bc, G, 'EAST', pz, 'sat', [1 0]); % 10 barsa

tot_time = 800*day(); % 200 days
n_steps = 50;
dt = repmat(tot_time/n_steps, [n_steps, 1]);

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

%% Initial condition
state = initResSol(G, 0*barsa, [1,0]);
t = 0;
p_mean = zeros(max(kk), n_steps);

f2 = figure(2); % to hold saturations

plotFaces(G, imperm_faces(:), 'FaceColor', 'none', 'EdgeColor', 'red', 'Linewidth', 1.5);
plotCellData(G, state.s(:,1), 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
colormap(flipud(winter)); colorbar; caxis([0, 1]);
title({'Saturation (1 -> water, 0 -> CO2)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f2, 'summer_sintef/case1/plots/sat_0', 'png');

f3 = figure(3); % to plot pressure

plotFaces(G, imperm_faces(:), 'FaceColor', 'none', 'EdgeColor', 'red', 'Linewidth', 1.5);
plotCellData(G, state.pressure, 'EdgeColor', 'none');
plotGrid(G, W.cells, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
colorbar;
title({'Pressure distribution (Pascal)' ['Time: ', formatTimeRange(t)]});
axis equal tight
view([0, 0])
drawnow

saveas(f3, 'summer_sintef/case1/plots/pres_0', 'png');

%% Copmute solutions
[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

dt_plot = 5; % plot every 5th iteration

for i=1:numel(states)
    t = t+dt(i);
    S = states{i}.s(:,1);
    assert(max(S) < 1+eps && min(S) > -eps);
    
    p = reshape(states{i}.pressure, [max(ii), max(kk)]);
    p_mean(:,i) = mean(p, 1).';
    
    figure(2);
    plotCellData(G, S, 'EdgeColor', 'none');
    colormap(flipud(winter)); colorbar; caxis([0 1]);
    title({'Saturation (1 -> water, 0 -> CO2)' ['Time:', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow, pause(0.2)
    
    figure(3);
    plotCellData(G, states{i}.pressure, 'EdgeColor', 'none');
    colorbar;
    title({'Pressure distribution (Pascal)' ['Time: ', formatTimeRange(t)]});
    axis equal tight
    view([0, 0])
    drawnow, pause(0.2)
    
    if ~mod(i, dt_plot)
       filename_f2 = sprintf('summer_sintef/case1/plots/sat_%d', i);
       filename_f3 = sprintf('summer_sintef/case1/plots/pres_%d', i);
       saveas(f2, filename_f2, 'png');
       saveas(f3, filename_f3, 'png');
   end
end

f4 = figure(4);
pcolor(p_mean); shading interp;
set(gca, 'YDir', 'reverse');
title('Horizontally averaged pressure');
xlabel('Time steps');
ylabel('Depth (m)');
colorbar;
drawnow;
%saveas(f4, 'summer_sintef/case1/plots/avg_pressure', 'png');

%% TESTING
G.cells.indexMap(k == 2);