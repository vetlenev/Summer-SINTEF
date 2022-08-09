classdef Hysteresis
    %HYSTERESIS To include hysteresis effects for CO2 drainage and water
    %imbibition
    
    methods
        function [krn] = Hysteresis(Sn, Sm1, Sm2, swr, snr, lambda)
          % Compute next relative permeability curves given current
          % drainage + imbibition curves.
          % Killough model used.
          
          % Input:
          %     Sn: CO2 saturation at current state
          %     Sm1: saturation at previous state
          %     Sm2: saturation at second previous state
          %     swr: residual water sat (endpoint of primary d-curve?)
          %     snr: residual CO2 sat (endpoint of primary i-curve?=
          %     lambda: interpolation coefficient
           
          % Output:
          %     krn: new relperm curve (for each cell in grid)
          Sw = 1 - Sn;
          if any(Sw > 1-snr) % only possible for at first displacement => we are at primary drainage
              Sw_scaled = (Sw - swr)/(1 - swr);
          else
              Sw_scaled = (Sw - swr)/(1 - snr - swr); % all subsequent drainage/imbibition
          end
          
          Sn_scaled = 1 - Sw_scaled;
          
          krn_D = primaryDrainage(Sn_scaled, swr, snr); % primary drainage curve
          krn_I = primaryImbibition(Sn_scaled, swr, snr); % primary imbibition curve
          S_D = inverseFunc(krn_D);
          S_I = inverseFunc(krn_I);
          S_nmax = max(S_D); % max from primary drainage (or just 1 - swr !??)
          S_nrmax = max(S_I); % max from primary imbibition ( or just snr !??)
          
          C = 1/S_nmax - 1/S_nrmax; % or 1/(1-swr) - 1/snr
          % Calculate S_ni (saturation where we go from drainage to
          % imbibition)
          drainage = ones(numel(Sn), 1); % boolean: 1 if drainage, 0 if imbibition
          drainage(Sn < Sn_m1 & Sn_m1 > Sn_m2) = 0; % CO2 saturation goes from increasing to decreasing in this cell -> imbibition starts!
          S_ni(ismember(drainage, 0)) = S(ismember(drainage, 0));
          % for the other elements, S_ni not known yet, and we still follow
          % primary drainage curve
          S_nr = 1 / (C + 1/S_ni);
          krn = krn_D(S_ni)*((S - S_nr)/(S_ni - S_nr)).^lambda; % new imbibition curve (equals subsequent drainage curve)          
      end          
          
    end
end

