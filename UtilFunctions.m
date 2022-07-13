classdef UtilFunctions
   methods (Static)
       function f = fullsizeFig(k)          
           f = figure(k);
           f.Position = get(0, 'Screensize');           
       end
       
       function sorted_struct = sortStructByField(my_struct, fieldname)
           tab = struct2table(my_struct);
           sorted_tab = sortrows(tab, fieldname);
           sorted_struct = table2struct(sorted_tab);
       end
       
      function pc_val = Pc(S, swr, p_e, cap, n)
         S_scaled = max(S-swr, 1e-5);
         pc_val = p_e*S_scaled.^(-1/n); % Corey model
         pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         pc_val(S<=swr) = cap;
         pc_val(S==1) = 0; % No pressure if no saturation
      end
      
      function pc_val = LeverettJ(S, swr, snr, phi, K, K_base, n, J1, J2, J3, k1, k2)                 
         S_scaled = (S - swr) / (1 - snr- swr);                 
         J = J1./((1 + k1*S_scaled).^n) - J2./((1 + k2*(1-S_scaled)).^n) + J3;  
         theta = 0;
         sigma = 0.02;            
         pc_val = sigma*cos(theta)*sqrt(phi./K).*J;           
         %pc_val(S <= swr) = max(pc_val(S > swr)); % cap value at residual saturation 
         %pc_val(S == 0) = 0; % No pressure if only water apparent
         pc_val(K == K_base) = 0; % No capillary pressure in background
      end
      
      function [tot_vol, region_vol] = Co2VolumeRatio(G, cells, S, rock, fluid)
          % Estimates volume of CO2 stored in provided cells.
          % Inputs:
          %     G: grid object
          %     cells: selected cells to estimate volume for (i.e. top)
          %     S: CO2 saturation values for each cell in grid at current state
          %     rock: rock object
          %     fluid: fluid object
          % Outputs:
          %     tot_vol: total volume of CO2 in grid
          %     region_vol: volume of CO2 in region represented by
          %     selected cells                   
          tot_vol = sum(G.cells.volumes.*S(G.cells.indexMap))*mean(rock.poro);%*fluid.rhoOS;               
          region_vol = sum(G.cells.volumes(cells).*S(cells))*mean(rock.poro);%*fluid.rhoOS;      
      end
      
      function [res_vol] = Co2Residual(G, S, sor, buff, S_gt_sor, rock)
          % Estimates volume of CO2 residually trapped in domain.
          % Inputs:
          %     G: grid object         
          %     S: CO2 saturation values for each cell in grid at current state        
          %     sor: residual oil saturation
          %     buff: small buffer for residual saturation
          %     S_gt_sor: boolean array for each cell, true if oil
          %     saturation has been larger than sor at any time step         
          %     rock: rock object
          %     fluid: fluid object
          % Outputs:
          %     res_vol: residual volume of CO2         
          S(S > sor+buff) = 0.0; % will take ages for oil sat to actually reach sor
          not_gt_sor = ~any(S_gt_sor);
          S(not_gt_sor(:) & S <= sor) = 0.0; % first occurence of sor is from injection, not residual saturation
          res_vol = sum(G.cells.volumes.*S)*mean(rock.poro);%*fluid.rhoOS;
      end
      
      function [res_vol] = Co2ResidualTrapped(G, S, sor, rock)
          % Estimates volume of CO2 residually trapped in domain.
          % Inputs:
          %     G: grid object         
          %     S: CO2 saturation values for each cell in grid at current state        
          %     sor: residual oil saturation        
          %     rock: rock object         
          % Outputs:
          %     res_vol: residual volume of CO2
          S_res = min(S, sor);         
          res_vol = sum(G.cells.volumes.*S_res)*mean(rock.poro);%*fluid.rhoOS;
      end
      
      function [struct_vol] = Co2StructuralTrapped(G, S, sor, trapped_cells, rock)
          % Estimates volume of CO2 structurally trapped in domain, defined
          % as all cells residing below concave shape but above bottommost
          % part of layer, including cells immediately below layer.
          % 
          % Inputs:
          %     G: grid object         
          %     S: CO2 saturation values for each cell in grid at current state        
          %     sor: residual oil saturation   
          %     trapped_cells: cell candidates for structral trapping
          %     rock: rock object         
          % Outputs:
          %     res_vol: residual volume of CO2
          
          %S_nonres = S(S > sor);
          non_trapped_cells = G.cells.indexMap(setdiff(G.cells.indexMap, trapped_cells));
          S(non_trapped_cells) = 0; % omit cells not candidated for structural trapping
          S_struct = max(S-sor, 0); % omit residually trapped CO2
          struct_vol = sum(G.cells.volumes.*S_struct)*mean(rock.poro);                       
      end
      
      function [S] = Hysteresis(Smax, swr, snr)
          % Compute next relative permeability curves given current
          % drainage + imbibition curves and maximum Smax (max saturation of
          % non-wetting phase) from history of states.
          krn = computeDrainageImbibitionCurves(); % for next state
          S = inverseFunc(krn);    
          Smax = max(max(S), Smax);
      end          
          
   end
end