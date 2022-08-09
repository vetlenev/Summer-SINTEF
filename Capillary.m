classdef Capillary
    %CAPILLARY Functions for capillary pressure
    
    methods (Static)
        function pc_val = Pc(S, swr, p_e, cap, n)
         S_scaled = max(S-swr, 1e-5);
         pc_val = p_e*S_scaled.^(-1/n); % Corey model
         pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         pc_val(S<=swr) = cap;
         pc_val(S==1) = 0; % No pressure if no saturation
      end
      
      function pc_val = PcNew(S, swr, snr, p_e, cap, n)
         S_scaled = max( (S-swr)/(1-snr-swr), 1e-5);
         pc_val = p_e*S_scaled.^(-1/n); % Corey model
         pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         pc_val(S<=swr) = cap;
         pc_val(S==1) = 0; % No pressure if no saturation
      end  
      
      function pc_val = LeverettJ(S, phi, K, K_base, median_pc)                 
         %S_scaled = max((S - swr) / (1 - snr- swr), 1e-5); % NB: no scaling of saturation for Leverett-J                                   
         surf_tension = median_pc/sqrt(median(phi)/median(K)); % median reservoir properties give cap pressure of median_pc (barsa)
         pc_val = surf_tension*sqrt(phi./K).*min(max((1 - S), 0), 1);         
         pc_val(K == K_base) = 0; % No capillary pressure in background
      end
      
      function pc_val = LeverettJnew(S, swr, snr, phi, K, K_base, n, J1, k1, k2)                 
         %S_scaled = (S - swr) / (1 - snr- swr);                         
         theta = 0;
         sigma = 0.02;                             
         J2 = J1 / (1 + k1);% to ensure J-func takes value zero at S = 1
         J = J1./((1 + k1*S_scaled).^n) - J2./((1 + k2*(1-S_scaled)).^n);% + J3; 
         pc_val = sigma*cos(theta)*sqrt(phi./K).*J; 
         pc_val(K == K_base) = 0; % No capillary pressure in background
      end
      
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
    end
end

