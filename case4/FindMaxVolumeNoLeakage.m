function [volumes, leaked_boundary, gt_leaked_perc, ...
            states, rates, scales, rel_diff, categorized_vols] = FindMaxVolumeNoLeakage(G, rock, rates, state, model, ...
                                                                                trapped_imperm, trapped_lowperm, sor, swr, ...
                                                                                dt, bc, inj_stop_ratio, leaked_perc, kmax)
% Run simulations for different injection rates, and return max injection
% rate that does not yield any leakage.
% Return total volume that can be injected without leakage, and what open
% boundary (top or right) caused the leakage
        
n_steps = numel(dt);
x = G.cells.centroids(:,1);
z = G.cells.centroids(:,3);
right_cells = G.cells.indexMap(x == max(x));
top_cells = G.cells.indexMap(z == min(z));
[ii, jj, kk] = gridLogicalIndices(G);
%right_faces = gridCellFaces(G, right_cells);
%top_faces = gridCellFaces(G, top_cells);
%right_faces = right_faces(2:6:numel(right_faces));
%top_faces = top_faces(4:6:numel(top_faces));
dims = G.cartDims;
nz = dims(3);
cc_cf_x = max(diff(x))/2; % distance from cell centroids to face centroids in x-direction
cc_cf_z = max(diff(z))/2; % NB: assumes uniform spacing!
lx = max(x) + cc_cf_x;
lz = max(z) + cc_cf_z;

volumes = [];
gt_leaked_perc = [1]; % boolean storing if leakage has occured for given injection rate
scales = []; scales_gt = []; scales_lt = [];
leaked_boundary = 'none';

k = 1;
epsilon = 0.01;
rel_diff = [2*epsilon];

categorized_vols = {{}, {}, {}, {}}; % residual, structural, free plume, leaked

%while (k <= kmax || gt_leaked_perc(end) == 1) && ( rel_diff(end) > epsilon || gt_leaked_perc(end) == 1 )
while k <= kmax && ( rel_diff(end) > epsilon || gt_leaked_perc(end) == 1 )
    if k == 1
        gt_leaked_perc = [];
    end
    t = 0;
    rate = rates(k);  
    
    well_h = 1; % cell perforations in vertical direction
    perforation_idx = G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
    W = addWell([], G, rock, perforation_idx, ...
            'Type', 'rate', 'Val', rate, ...
            'InnerProduct', 'ip_tpf', ...
            'Radius', 0.1, 'Dir', 'x', ...
            'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
            'Name', 'P1');  
        
    schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
    times = cumsum(dt)/year();
    %inj_stop = fix(inj_stop_ratio*n_steps);
    [~, inj_stop] = min(abs(times - inj_stop_ratio*times(end)));
    schedule.control(2) = schedule.control(1); % create second well
    schedule.control(2).W.status = 0; % shut off second well
    schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway      
    
    [wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);
        
    tot_vols = zeros(numel(states)+1, 1);
    simulation_vols = zeros(size(tot_vols));
    
    residual_vols = zeros(size(tot_vols)); % first element is initial state    
    structural_vols = zeros(size(tot_vols));
    structural_res_vols = zeros(size(tot_vols));
    free_vols = zeros(size(tot_vols));
    leaked_vols = zeros(size(tot_vols));   
    leaked_ratios = zeros(size(tot_vols));
    
    for i=1:numel(states)       
        t = t + dt(i);
        Sw = states{i}.s(:,1);
        Sn = 1 - Sw;
        assert(max(Sw) < 1+eps && min(Sw) > -eps);               
               
        simulation_vols(i+1) = min(rate*t, rate*sum(dt(1:inj_stop-1)));       
        tot_vols(i+1) = min(sum(G.cells.volumes.*Sn)*mean(rock.poro), simulation_vols(i+1)); % prevent cases where tot_vol is larger than sim_vol (impossible due to no source terms)
        
        if i >= inj_stop-1
            residual_vols(i+1) = VolumeTrapping.Co2ResidualTrapped(G, Sn, sor, rock);
        end
        % only include volume trapped under imperm layers
        [structural_vols(i+1), ~, ...
          structural_res_vols(i+1), ~] = VolumeTrapping.Co2StructuralTrapped(G, Sn, sor, trapped_imperm, trapped_lowperm, rock); 
        leaked_vols(i+1) = simulation_vols(i+1) - tot_vols(i+1);
        leaked_ratios(i+1) = leaked_vols(i+1) / simulation_vols(i+1);      
    end    
    
    disp(min(leaked_ratios))
    disp(max(leaked_ratios))
    
    % linear interpolation of residual trapping (assume gradual increase in
    % imbibition over time until injection stop
    dt_samples = [1, inj_stop];
    res_samples = residual_vols(dt_samples);
    dt_interp = 1:inj_stop;
    residual_vols(1:inj_stop) = interp1(dt_samples, res_samples, dt_interp, 'linear');
    
    for i=1:numel(states) % Need to caclulate free plume AFTER interpolation of residual
        free_vols(i+1) = max(0, tot_vols(i+1) - (residual_vols(i+1) + structural_vols(i+1))); % interpolation of residual vol could cause negative value at start --> prevent this!
    end
    
    volumes = cat(1, volumes, simulation_vols(end)); % store total volume from previous simulation
    categorized_vols{1} = residual_vols;
    categorized_vols{2} = structural_vols;
    categorized_vols{3} = free_vols;
    categorized_vols{4} = leaked_vols;
    
    % Calculate leakage and compute new rates
    %leakage = any(leaked_ratios > leaked_perc + 1e-6); % avoid round-off errors / negligable "particle" leakage
    % OR: Only need to check leaked rate at end of simulation, as this is
    % guaranteed to be the maximum leaked vol achieved during simulation
    leakage = leaked_ratios(end) > leaked_perc + 1e-6;   
    
    if leakage
        if any(Sn(right_cells)) && any(Sn(top_cells))
            leaked_boundary = ['right' 'top'];
        elseif any(Sn(right_cells))
            leaked_boundary = ['right'];          
        elseif any(Sn(top_cells))
            leaked_boundary = ['top'];    
        end        
    end
       
    gt_leaked_perc = cat(1, gt_leaked_perc, leakage); % boolean indicating if simulation yielded leakage greater than provided percentage 
    idx_last_rates = {find(ismember(gt_leaked_perc, 0)), ...
                        find(ismember(gt_leaked_perc, 1))};    
    num_lt = numel(idx_last_rates{1}); % number of simulations with leakage LESS than perc
    num_gt = numel(idx_last_rates{2}); % num simulations with leakage GREATER than perc
    
    % set new rate based on how early leakage occured
    %scale_gt = leaked_ratios(end);
    scale_gt = numel(find(leaked_ratios > leaked_perc + 1e-5)) / numel(leaked_ratios);
    scales_gt = cat(1, scales_gt, scale_gt);
    
    % set new rate based on how close CO2 is to reach leakage perc
    Sn_closest_right = max(ii(Sn > 1e-6));
    Sn_closest_top = min(kk(Sn > 1e-6));
    % select minimum, as this is the side closest to leakage
    scale_lt = min( (max(ii) - Sn_closest_right)/max(ii), ... 
                    Sn_closest_top/max(kk) ) + leaked_perc*exp(1+leaked_perc); % last term adds more senitivity to large leakage percentages, because we need larger injection rates for larger leakage threshold
    scales_lt = cat(1, scales_lt, scale_lt);      
                
    if (num_lt && num_gt) % choose middle point as new rate
        idx_last_lt = idx_last_rates{1}(end); % index of last rate BELOW leakage perc
        idx_last_gt = idx_last_rates{2}(end); % index of last rate ABOVE leakage perc
        new_rate = 0.5*(rates(idx_last_lt) + rates(idx_last_gt));
       
    elseif num_lt % lower than perc -> exponentially increase previous rate
        idx_last_lt = idx_last_rates{1}(end);
        disp(scale_lt)        
        new_rate = exp(scale_lt*k)*rates(k);
    else % above perc -> exponentially decrease previous rate
        %idx_last_gt = idx_last_rates{2}(end);
        %new_rate = (1-scale_gt)*rates(k)*exp(-scale_gt*k);    
        new_rate = rates(k)*exp(-scale_gt*k);
    end
      
    rates = cat(1, rates, new_rate); % NB: unit m^3/s
    
    rel_diff = cat(1, rel_diff, abs(rates(end) - rates(end-1))/rates(end));
    disp('Done with iteration nr:')
    disp(k);
    disp('New rate:')
    disp(new_rate)
    disp('Leaked:')
    disp(leakage)
    k = k+1;   
end

rel_diff(1) = []; % delete first dummy value
rates(end) = []; % last rate not used

% Compute utilized percentage of (permanently) structural traps (taking
% into account residual saturation!)
structural_utilized = (structural_vols(end) + structural_res_vols(end)) / (sum(G.cells.volumes(trapped_imperm))*mean(rock.poro)*(1-swr));
disp('Percentage of structural traps utilized')
disp(structural_utilized*100)

if k == numel(rates) && strcmp(leaked_boundary, 'none')
    disp('No leakage for selected injection rates. Consider increasing the rates to find even better solution');
elseif k > kmax && (gt_leaked_perc(end) == 1)
    disp('Max iteration reached, and leakage still exceeds the prescribed percentage.');
end

end

