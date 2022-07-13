function [volumes, leaked_boundary, gt_leaked_perc, ...
            states, rates, scales, rel_diff, categorized_vols] = FindMaxVolumeNoLeakage(G, rock, rates, state, model, ...
                                                                                trapped_cells, sor, dt, bc, ...
                                                                                inj_stop_ratio, leaked_perc, kmax)
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

categorized_vols = {{}, {}, {}, {}}; % residual, structural, free plume

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
    
    inj_stop = fix(inj_stop_ratio*n_steps);
    schedule.control(2) = schedule.control(1); % create second well
    schedule.control(2).W.status = 0; % shut off second well
    schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway      
    
    [wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);
    leaked_ratios = zeros(numel(states)+1, 1);
    residual_vols = zeros(numel(states)+1, 1); % first element is initial state
    structural_vols = zeros(numel(states)+1, 1);
    free_vols = zeros(numel(states)+1, 1);
    simulation_vols = zeros(numel(states)+1, 1);
    
    for i=1:numel(states)       
        t = t + dt(i);
        Sw = states{i}.s(:,1);
        Sn = 1 - Sw;
        assert(max(Sw) < 1+eps && min(Sw) > -eps);               
        
        tot_vol = sum(G.cells.volumes.*Sn)*mean(rock.poro);
        simulation_vols(i+1) = min(rate*t, rate*sum(dt(1:inj_stop-1)));
        
        if i >= inj_stop-1
            residual_vols(i+1) = UtilFunctions.Co2ResidualTrapped(G, Sn, sor, rock);
        end
        structural_vols(i+1) = UtilFunctions.Co2StructuralTrapped(G, Sn, sor, trapped_cells, rock);        
        leaked_ratios(i+1) = (simulation_vols(i+1) - tot_vol) / simulation_vols(i+1);
    end
    
    % linear interpolation of residual trapping
    dt_samples = [1, inj_stop];
    res_samples = residual_vols(dt_samples);
    dt_interp = 1:inj_stop;
    residual_vols(1:inj_stop) = interp1(dt_samples, res_samples, dt_interp, 'linear');
    
    for i=1:numel(states) % Need to caclulate free plume AFTER interpolation of residual
        free_vols(i+1) = max(0, simulation_vols(i+1) - (residual_vols(i+1) + structural_vols(i+1))); % interpolation of residual vol could cause negative value at start --> prevent this!
    end
    
    volumes = cat(1, volumes, simulation_vols(end)); % store total volume from previous simulation
    categorized_vols{1} = residual_vols;
    categorized_vols{2} = structural_vols;
    categorized_vols{3} = free_vols;
    categorized_vols{4} = leaked_ratios;
    
    % Calculate leakage and compute new rates
    leakage = any(leaked_ratios > leaked_perc + 1e-5); % avoid round-off errors / negligable "particle" leakage  
    
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
    Sn_closest_right = max(ii(Sn > 1e-5));
    Sn_closest_top = min(kk(Sn > 1e-5));
    % select minimum, as this is the side closest to leakage
    scale_lt = min( (max(ii) - Sn_closest_right)/max(ii), ... 
                    Sn_closest_top/max(kk) ) + leaked_perc*exp(1+leaked_perc); % last term adds more senitivity to large leakage percentages
    scales_lt = cat(1, scales_lt, scale_lt);
                
    scales = cat(1, scales, [scale_lt, scale_gt]);
    disp('scales')
    disp(scales)
                
    if (num_lt && num_gt) % choose middle point as new rate
        idx_last_lt = idx_last_rates{1}(end); % index of last rate BELOW leakage perc
        idx_last_gt = idx_last_rates{2}(end); % index of last rate ABOVE leakage perc
        new_rate = 0.5*(rates(idx_last_lt) + rates(idx_last_gt));
%         rate_lt = rates(idx_last_lt) / (rates(idx_last_lt) + rates(idx_last_gt)); % x_A in book
%         rate_gt = rates(idx_last_gt) / (rates(idx_last_lt) + rates(idx_last_gt)); % x_B in book
%         mid_shift = (scales_lt(idx_last_lt) - scales_gt(idx_last_gt))/(scales_lt(idx_last_lt) ...
%                                                     + scales_gt(idx_last_gt)) * 0.5*(rate_gt - rate_lt);
%         new_scale = 0.5*(rate_lt + rate_gt) + mid_shift;
%         new_rate = new_scale * (rates(idx_last_lt) + rates(idx_last_gt));        
    elseif num_lt % lower than perc -> exponentially increase previous rate
        idx_last_lt = idx_last_rates{1}(end);
        disp(idx_last_lt)        
        new_rate = exp(scale_lt*k)*rates(k);
    else % above perc -> exponentially decrease previous rate
        %idx_last_gt = idx_last_rates{2}(end);
        new_rate = -exp(scale_gt*k)*rates(k);    
    end
      
    rates = cat(1, rates, new_rate);

    %disp(rates)
    rel_diff = cat(1, rel_diff, abs(rates(end) - rates(end-1))/rates(end));
    disp('Done with iteration nr:')
    disp(k);
    disp('New rate:')
    disp(new_rate)
    k = k+1;   
end

rel_diff(1) = []; % delete first dummy value

if k == numel(rates) && strcmp(leaked_boundary, 'none')
    disp('No leakage for selected injection rates. Consider increasing the rates to find even better solution');
elseif k > kmax && (gt_leaked_perc(end) == 1)
    disp('Max iteration reached, and leakage still exceeds the prescribed percentage.');
end

end

