classdef GridSetupParams < handle
    %GRIDSETUP Functions for generating grid, tailored for sensitivity
    %analysis for parameters (e.g. correlation length).
    %parameters 
    properties
        G
        X
        Z
        nLowperm
        nImperm
        ii
        kk
        dX
        angle
        xStart
        xStop
        zStart
        zStop
        log_min
        log_max
        corr_len_x
        corr_len_z
        Perm
        lowperm
        baseperm
    end
    
    methods
        function o = GridSetupParams(G, n_lowperm_layers, n_imperm_layers, perms)
            [ii, ~, kk] = gridLogicalIndices(G);            
            
            if nargin == 4
                o.G = G;
                o.X = G.cells.centroids(:,1);
                o.Z = G.cells.centroids(:,3);
                o.nLowperm = n_lowperm_layers;
                o.nImperm = n_imperm_layers;
                o.ii = ii;
                o.kk = kk;
                o.dX = mode(diff(o.X));
                
                o.angle = pi/4;
                o.xStart = struct; o.xStop = struct;
                o.zStart = struct; o.zStop = struct;               
                o.corr_len_x = 1; o.corr_len_z = 1;
                o.log_min = 0.01; o.log_max = 100;
                           
                o.lowperm = perms{1};
                o.baseperm = perms{2};
                o.Perm = repmat(o.baseperm, [o.G.cells.num 1]);
            end
        end
        
        function [line_idx, anticline_idx] = DefineLayers(o)
            % Randomly distribute bunch of low/imperm layers           
            
            perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);         

            anticline_idx = 1:3:(o.nLowperm+o.nImperm); % every third layer is anticline
            line_idx = setdiff(1:(o.nLowperm+o.nImperm), anticline_idx); % indices for line layers
            
            rand_x_mid = linspace(min([o.X])-o.dX/2, max(o.X)+o.dX/2, o.nLowperm+o.nImperm); % o.dX/2 to move from centroid to endpoint
            rand_x_mid = rand_x_mid(randperm(length(rand_x_mid))); % random permute

            o.xStart.lowperm = rand_x_mid(1:o.nLowperm);
            o.xStart.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm);

            x_length_lowperm = randi([round(perc_val(o.X, 0.05)), ...
                                                round(perc_val(o.X, 0.35))], ...
                                                [o.nLowperm, 1]).';

            o.xStart.lowperm = rand_x_mid(1:o.nLowperm) - round(x_length_lowperm/2);
            o.xStop.lowperm = rand_x_mid(1:o.nLowperm) + round(x_length_lowperm/2);

            x_length_imperm = randi([round(perc_val(o.X, 0.05)), ...
                                                round(perc_val(o.X, 0.35))], ...
                                                [o.nImperm, 1]).';

            o.xStart.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm) - round(x_length_imperm/2);
            o.xStop.imperm = rand_x_mid(o.nLowperm+1:o.nLowperm+o.nImperm) + round(x_length_imperm/2);           

            rand_z_start = linspace(perc_val(o.Z, 0.05), perc_val(o.Z, 0.9), o.nLowperm+o.nImperm);
            o.zStart.lowperm = sort(randsample(rand_z_start, o.nLowperm));
            rand_z_start(ismember(rand_z_start, o.zStart.lowperm)) = [];
            o.zStart.imperm = sort(rand_z_start);

            o.zStop.lowperm = o.zStart.lowperm + perc_val(o.Z, 0.02);
            o.zStop.imperm = o.zStart.imperm + perc_val(o.Z, 0.02); 
            
            o.corr_len_x = mean(o.xStop.lowperm - o.xStart.lowperm) / 100; % default value
            o.corr_len_z = mean(o.zStop.lowperm - o.zStart.lowperm) / 10; % default value
            
        end
        
        function [added_layers, trapped_cells, ...
                    line_idxs, anticline_idxs, ...
                    all_anticline_dummy] = ...
                GenerateLayers(o, line_idx, anticline_idx) 
            % Include defined layers in grid and set their physical
            % properties
            perc_val = @(arr, perc, varargin) perc*(max(arr) - min(arr)) + min(arr);
                        
            layer_types = {'imperm', 'lowperm'};

            anticline_idx_lowperm = anticline_idx(1:round(o.nLowperm/(o.nLowperm+o.nImperm)*length(anticline_idx)));
            anticline_idx_lowperm(anticline_idx_lowperm > o.nLowperm) = []; % remove invalid values
            anticline_idx_imperm = anticline_idx(find(anticline_idx == anticline_idx_lowperm(end))+1:end) - o.nLowperm;
            anticline_idx_imperm(anticline_idx_imperm > o.nImperm) = []; % remove invalid values

            anticline_idxs = {anticline_idx_imperm, anticline_idx_lowperm};

            line_idx_lowperm = line_idx(1:round(o.nLowperm/(o.nLowperm+o.nImperm)*length(line_idx)));
            line_idx_lowperm(line_idx_lowperm > o.nLowperm) = []; % remove invalid values
            line_idx_imperm = line_idx(find(line_idx == line_idx_lowperm(end))+1:end) - o.nLowperm;
            line_idx_imperm(line_idx_imperm > o.nImperm) = []; % remove invalid values

            line_idxs = {line_idx_imperm, line_idx_lowperm};            
            
            o.xStart.imperm(anticline_idxs{1}) = o.xStart.imperm(anticline_idxs{1}) - perc_val(o.X, 0.05);
            o.xStop.imperm(anticline_idxs{1}) = o.xStop.imperm(anticline_idxs{1}) + perc_val(o.X, 0.05);
            o.xStart.lowperm(anticline_idxs{2}) = o.xStart.lowperm(anticline_idxs{2}) - perc_val(o.X, 0.05);
            o.xStop.lowperm(anticline_idxs{2}) = o.xStop.lowperm(anticline_idxs{2}) + perc_val(o.X, 0.05);           
            
            o.zStop.imperm(anticline_idxs{1}) = o.zStop.imperm(anticline_idxs{1}) + perc_val(o.Z, 0.01);
            o.zStop.lowperm(anticline_idxs{2}) = o.zStop.lowperm(anticline_idxs{2}) + perc_val(o.Z, 0.01);                                                                     
            
            added_layers = {cell(1,numel(line_idxs{1})+numel(anticline_idxs{1})), ...
                            cell(1,numel(line_idxs{2})+numel(anticline_idxs{2}))};
                        
            all_anticline_dummy = {cell(1, numel(anticline_idxs{1})), ...
                                    cell(1, numel(anticline_idxs{2}))};          
            
            trapped_cells = {{}, {}}; % {structural imperm, structural free}

            for k=1:numel(layer_types)
                layer_type = layer_types{k};
                idx_count = 1;
                % Make straight forms
                for i=line_idxs{k}
                    if i == line_idxs{k}(fix(end/2)) && strcmp(layer_type, 'lowperm') % Ensure one far-extending layer
                        o.xStart.(layer_type)(i) = min(o.X);
                        o.xStop.(layer_type)(i) = max(o.X);
                        o.zStop.(layer_type)(i) = o.zStart.(layer_type)(i) + perc_val(o.Z, 0.04); % increase vertical extent of the full-blocking layer                        
                        
                        line_cells = o.G.cells.indexMap(o.X >= o.xStart.(layer_type)(i) & ...
                                                    o.X <= o.xStop.(layer_type)(i) & ...
                                                    o.Z > o.zStart.(layer_type)(i) & ...
                                                    o.Z < o.zStop.(layer_type)(i));
                                                
                        [hole_min, hole_max] = deal(0.05, 0.15); % min/max length of vertical extent
                        line_cells = UtilFunctions.cutLayer(o, line_cells, hole_min, hole_max, layer_type, i);
                        
                    else
                        line_cells = o.G.cells.indexMap(o.X >= o.xStart.(layer_type)(i) & ...
                                                    o.X <= o.xStop.(layer_type)(i) & ...
                                                    o.Z > o.zStart.(layer_type)(i) & ...
                                                    o.Z < o.zStop.(layer_type)(i)); 
                                                
                    end
                    
                    added_layers{k}{idx_count} = line_cells;
                    idx_count = idx_count + 1;
                end

                % Make anticline forms
                for i=anticline_idxs{k}
                    num_x_cells = numel(o.G.cells.indexMap(o.X > o.xStart.(layer_type)(i) ... 
                                        & o.X < o.xStop.(layer_type)(i) ...
                                        & o.Z == min(o.Z))); % only select ONE horizontal patch
                    theta0 = o.angle;
                    theta = linspace(theta0, pi - theta0, num_x_cells);
                    r = (o.xStop.(layer_type)(i) - o.xStart.(layer_type)(i))/2;
                    z_anticline_start = -r*sin(theta) + o.zStart.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % -r*sin to get anticline (since z positive downwards)
                    z_anticline_stop = -r*sin(theta) + o.zStop.(layer_type)(i) + r/2*(1+cos(pi/2-theta0)); % + r*cos(pi/4) to shift to bottom of spherical cap

                    anticline_cells = [];
                    anticline_cells_dummy = o.G.cells.indexMap(o.X > o.xStart.(layer_type)(i) & ...
                                                    o.X < o.xStop.(layer_type)(i) & ...
                                                    o.Z > min(z_anticline_start) & ...
                                                    o.Z < max(z_anticline_stop));

                    max_z = max(o.kk(anticline_cells_dummy));
                    min_z = min(o.kk(anticline_cells_dummy)); 
                    z_bottom = unique(o.Z(o.kk == max_z)); % z-coord bottom part of layer                    

                    trapped_cells_cat = [];

                    for j=1:num_x_cells
                        x_slice_start = (o.xStop.(layer_type)(i) - o.xStart.(layer_type)(i))*(j-1)/num_x_cells + o.xStart.(layer_type)(i);
                        x_slice_stop = (o.xStop.(layer_type)(i) - o.xStart.(layer_type)(i))*j/num_x_cells + o.xStart.(layer_type)(i);
                        anticline_single = o.G.cells.indexMap(o.X > x_slice_start & ...
                                                    o.X < x_slice_stop & ...
                                                    o.Z > z_anticline_start(j) & ...
                                                    o.Z < z_anticline_stop(j));
                                                
                        trapped_slice = o.G.cells.indexMap(o.X > x_slice_start & ...
                                                     o.X < x_slice_stop & ...
                                                     o.Z > z_anticline_stop(j) & ...
                                                     o.Z <= z_bottom).';
                                                 
                        trapped_cells_cat = cat(2, trapped_cells_cat, trapped_slice);

                        anticline_cells = cat(1, anticline_cells, anticline_single);
                    end

                    trapped_cells{k} = cat(1, trapped_cells{k}, trapped_cells_cat);
                    
                    added_layers{k}{idx_count} = anticline_cells;
                    all_anticline_dummy{k}{idx_count-numel(line_idxs{k})} = anticline_cells_dummy;
                    idx_count = idx_count + 1;
                end
                
                trapped_cells{k} = horzcat(trapped_cells{k}{:});
            end
                        
        end
        
        function [trapped_cells] = SetPermeability(o, added_layers, trapped_cells, ...
                                                        line_idxs, anticline_idxs, ...
                                                        anticline_cells_dummy)
            perm = repmat(o.baseperm, [o.G.cells.num 1]);
            perm_vals = {1e-5*milli*darcy, o.lowperm};
            
            layer_types = {'imperm', 'lowperm'};
            new_added_layers = {[], []};        
            
            for k=1:numel(layer_types)
                layer_type = layer_types{k};               
                % Make straight forms
                for i=1:numel(line_idxs{k})
                    line_cells = added_layers{k}{i};
                    
                    [ix, ~, kz] = gridLogicalIndices(o.G, line_cells);
                    nxi = numel(unique(ix)); nzi = numel(unique(kz)); % dimensions of particular low-perm layer

                    line_cells_G = intersect(o.G.cells.indexMap, line_cells);
                    added_layers_G = intersect(o.G.cells.indexMap, vertcat(new_added_layers{:})); 
                    overlapped_cells = o.G.cells.indexMap(intersect(line_cells_G, ...
                                                                    added_layers_G));
                    store_perm = perm(overlapped_cells);
                    
                    %perm(line_cells) = perm_vals{k} + perm_vals{k}*FastGaussian([nxi nzi], 0.8, [o.corr_len_x o.corr_len_z])*strcmp(layer_type, 'lowperm');                    
                    log_field = exp(gaussianField([nxi nzi], ...
                                                   [log(o.log_min), log(o.log_max)], ...
                                                   [o.corr_len_x, 1, o.corr_len_z], 5));                   
                    log_field = reshape(log_field, nxi*nzi, 1);                   
                    %log_field = log_field - mean(log_field);
                  
                    if strcmp(layer_type, 'lowperm')
                        perm(line_cells) = perm_vals{k}*log_field;
                    else
                        perm(line_cells) = perm_vals{k};
                    end
                    perm(overlapped_cells) = store_perm;  
                    
                    new_added_layers{k} = cat(1, new_added_layers{k}, line_cells);
                end
                
                % Anticline forms
                for i=1:numel(anticline_idxs{k})
                    anticline_cells = added_layers{k}{numel(line_idxs{k})+i};
                    
                    [ix, ~, kz] = gridLogicalIndices(o.G, anticline_cells);
                    nxi = numel(unique(ix)); nzi = numel(unique(kz));
                                        
                    added_layers_G = intersect(intersect(o.G.cells.indexMap, anticline_cells_dummy{k}{i}), ...
                                                intersect(o.G.cells.indexMap, vertcat(new_added_layers{:}))); 
              
                    overlapped_cells = o.G.cells.indexMap(added_layers_G);
                    
                    store_perm = perm(overlapped_cells);                                        
                    
                    %log_field = exp(FastGaussian([nxi nzi], 0.8, [o.corr_len_x o.corr_len_z]));
                    %perm(anticline_cells_dummy{k}{i}) = perm_vals{k} + perm_vals{k}**strcmp(layer_type, 'lowperm'); % rectangular box around anticline
                    log_field = exp(gaussianField([nxi nzi], ...
                                                   [log(o.log_min), log(o.log_max)], ... % o.lowperm/o.baseperm
                                                   [o.corr_len_x, 1, o.corr_len_z], 5));                  
                    log_field = reshape(log_field, nxi*nzi, 1);
                    %log_field = log_field - mean(log_field);
                    if strcmp(layer_type, 'lowperm')
                        perm(anticline_cells_dummy{k}{i}) = perm_vals{k}*log_field;
                    else
                        perm(anticline_cells_dummy{k}{i}) = perm_vals{k};
                    end
                    
                    perm(overlapped_cells) = store_perm;
                    background_cells = o.G.cells.indexMap(setdiff(anticline_cells_dummy{k}{i}, ...
                                                        vertcat(anticline_cells, vertcat(new_added_layers{:}))));
                    perm(background_cells) = o.baseperm; % reset perm for parts of rectangle NOT part of anticline
                    
                    new_added_layers{k} = cat(1, new_added_layers{k}, anticline_cells);
                end
            end                       
            
            imperm_layers = new_added_layers{1};
            lowperm_layers = new_added_layers{2};
            
            clear new_added_layers;
                      
            % Include cells blocked by lowperm cells at left, right and top
            % in list of temporary trapped cells
            for c=1:numel(lowperm_layers)
                cell = lowperm_layers(c);
                
                c_neighbors = o.G.faces.neighbors;
                [idx_neighbors, ~] = find(c_neighbors == cell);
                c_neighbors = c_neighbors(idx_neighbors, :);
                
                blocking_neighbors = sort(c_neighbors(c_neighbors > 0 & c_neighbors ~= cell)); % sort to recognise neighbors (due to logical indexing of grid)               

                x_cell = unique(o.X(o.ii == o.ii(cell)));
                z_cell = unique(o.Z(o.kk == o.kk(cell)));
                open_boundary_cell = (x_cell == max(o.X) || z_cell == min(o.Z));
                
                background_field = setdiff(o.G.cells.indexMap, union(trapped_cells{2}, lowperm_layers));
                
                % check if blocking layers are all lowperm layers
                if (all(ismember(blocking_neighbors(1:end-1), lowperm_layers)) ... % don't include end -> this is bottom
                        || ~any(ismember(blocking_neighbors(1:end-1), background_field))) ...
                        && ~open_boundary_cell % cells at open boundary are not trapped
                    trapped_cells{2} = horzcat(trapped_cells{2}, cell);
                end
            end
            
            % Remove trapped cells overlapping with imperm layers
            for k=1:numel(layer_types)                
                overlap_trap_imperm = intersect(intersect(o.G.cells.indexMap, trapped_cells{k}), ...
                                                intersect(o.G.cells.indexMap, imperm_layers)); % added_layers{1} = imperm layers
                
                trapped_cells{k}(ismember(trapped_cells{k}, overlap_trap_imperm)) = [];
            end               

            %all_added_layers = vertcat(added_layers{:});            
            perm = perm.*sign(perm); % set all positive values
            perm = max(perm, perm_vals{1}); % perm must be positive -> cap at non-zero value to avoid singular matrix
            o.Perm = perm;
        end
    end
end

