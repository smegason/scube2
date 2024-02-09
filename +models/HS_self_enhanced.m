% simple "phenomenological" model that models only two species:
% Hedgehog and Scube
% Scube changes the local D/gamma and flux at the origin
% scube doesnt modulate the degradation rate bc that's complicated

classdef HS_self_enhanced < models.model
   properties
       a_H
       c_H1
       c_H2
       D_H
       a_S
       c_S
       K_H
       n_H
       D_S
       K_D
       D_HS_rel
       Phi_HS_rel
   end
   
   methods
       function o = HS_self_enhanced(paramsModel, paramsSys)
            % call superclass constructor
            o = o@models.model(paramsModel, paramsSys);
            
            o.numSpecies = 2;
            
            % unpack model parameters into class
            o.a_H = paramsModel.a_H;
            o.c_H1 = paramsModel.c_H1;
            o.c_H2 = paramsModel.c_H2;
            o.D_H = paramsModel.D_H;
            o.a_S = paramsModel.a_S;
            o.c_S = paramsModel.c_S;
            o.K_H = paramsModel.K_H;
            o.n_H = paramsModel.n_H;
            o.D_S = paramsModel.D_S;
            o.K_D = paramsModel.K_D;
            o.D_HS_rel = paramsModel.D_HS_rel;
            o.Phi_HS_rel = paramsModel.Phi_HS_rel;
        end

        function ydot = derivatives_0D(o, t, y)
            % local derivatives
            H = y(1, :);            
            S = y(2, :);
            
            dHdt = -o.c_H1 * H - o.c_H2 * H .^ 2;
            dSdt = o.a_S * (o.K_H ^ o.n_H) ./ (o.K_H ^ o.n_H + H .^ o.n_H) - o.c_S * S;
            
            ydot = [dHdt; dSdt];
        end
        
        function ydot = derivatives_diffusion(o, t, y)
            % spatial derivatives
            H = y(1, :);            
            S = y(2, :);
            
            % Hedgehog D and flux depend on scube concentration
            varphi = (S ./ (S + o.K_D));
            D = o.D_H * (1 + (o.D_HS_rel - 1) * varphi);
            Phi = o.a_H * (1 + (o.Phi_HS_rel - 1) * varphi);
            
            % scale diffusion constants appropriately
            dHdt = D / o.delta^2 .* o.laplacian(H, 'r', 'r');
            dSdt = o.D_S / o.delta^2 * o.laplacian(S, 'r', 'r');
            
            % flux
            dHdt(1) = dHdt(1) + Phi(1) / o.delta;
            
            ydot = [dHdt; dSdt];
        end
        
        function A = get_perception(o, y)
            H = y(1, :);
            
            % scube expression domain
            A = (o.K_H ^ o.n_H) ./ (o.K_H ^ o.n_H + H .^ o.n_H);
            
            % total Hedgehog
            % A = squeeze(H);
        end
        
        function initial_conditions = get_initial_conditions_0D(o)
            initial_conditions = [0, 0];
        end
        
        function value = get_dPrevCell(o)
            value = [1, 0; 0, 1];
        end
        
        function value = get_dNextCell(o)
            value = [1, 0; 0, 1];
        end
        
        function value = get_dCurrCell(o)
            value = [1, 1; 1, 1];
        end
    end       
end