% Model base class
% Sets up a set of differential equations to model the 
% time-evolution of a system with multiple diffusible
% species

% Things that should be specified
% derivatives_0D: "non-spatial" part of the DEs
% derivatives_diffusion: 'spatial" part of the DEs *
% get_perception: returns a row vector with perception
% get_initial_conditions_0D: initial conditions
% get_Ss_species_0D: indices of steady-state species during preequilibrium
% get_Ss_species_1D: indices of steady-state species
% get_dPrevCell: jacobian sparsity matrix between this cell and previous
% get_dCurrCell: jacobian sparsity matrix for each cell and itself
% get_dNextCell: jacobian sparsity matrix between this cell and next**

% *both these functions take a y_matrix - size (numSpecies, N) as input
% and output
% **in most cases (for diffusible species); a diagonal matrix, with a one
% if the species is diffusible

% Inputs to constructor:
%   paramsModel : structure specifying model parameters
%   paramsSys   : structure specifying "system" parameters 
% Constructor extracts parameters from parameter files
% get_jacobian_matrix computes sparsity pattern for Jacobian of 1D model
% get_Ss_species_1D specifies indicies of species to be tested for steady stateness for 1D model 
% get_Ss_species_0D specifies indices of species to be tested for steady stateness for model during pre-equilibration (0D)
% get_initial_conditions_0D specifies initial conditions before the pre-equilibration stage
% 
% model_0D or model_1D can be passed directly to an ODE solver
% model_0D and model_1D take t, y as input and returns ydot
%

classdef model
    properties (SetAccess=protected)
        N  % how many cells (spatial resolution)
    	numSpecies
        L
        delta  % size of one cell
    end
    
    methods
        function obj = model(paramsModel, paramsSys)
            % extract parameters from paramsModel and paramsSys struct
            % files
            
            obj.N = paramsSys.N;
            obj.L = paramsSys.L;
            
            % size of one cell
            obj.delta = obj.L / obj.N;
        end
        
        function ydot = laplacian(obj, y, leftBC, rightBC)
            % returns the result of applying the "operator"
            % d^2/dx^2 on the space-function y(x), subject to 
            % left and right BCs which can either be reflective or
            % absorptive
            
            ydot = zeros(size(y));
            
            switch nargin
                case 2  % use reflective BCs by default
                    leftBC = 'reflective';
                    rightBC = 'reflective';
                case 4
                    % do nothing
                otherwise
                    error('laplacian takes 1 or 3 arguments')
            end
            
            % implement BCs
            if strcmp(leftBC(1), 'r')
                ydot(1) = y(2) - y(1);
            else  % absorptive BC
                ydot(1) = y(2) - 2 * y(1);
            end
            
            if strcmp(rightBC(1), 'r')
                ydot(obj.N) = y(obj.N - 1) - y(obj.N);
            else  % absorptive BC
                ydot(obj.N) = y(obj.N - 1) - 2 * y(obj.N);
            end
            
            % interior regions; discretize laplacian
            ydot(2:obj.N-1) = y(3:obj.N) + y(1:obj.N-2) - 2 * y(2:obj.N-1);
        end

        function ydot = derivatives_0D(obj, y)
            % computes the "local" derivative of each species
            ydot = zeros(size(y));  % replace this in child classes
        end
        
        function ydot = derivatives_diffusion(obj, y)
            % computes spatial derivatives of each species
            ydot = zeros(size(y));  % replace this in child classes
        end
        
        function A = get_perception(obj, y)
            % given the y-matrix, returns the quantity to plot (perception)
            A = zeros(1, obj.N);
        end
        
        function ydot = model_0D(obj, t, y)
            % pass to ode solver
            ydot = derivatives_0D(t, y);
        end
        
        function ydot = model_1D(obj, t, y)
            % pass to ode solver
            y_matrix = reshape(y, obj.numSpecies, obj.N);
            ydot = obj.derivatives_0D(t, y_matrix) + obj.derivatives_diffusion(t, y_matrix);
            ydot = reshape(ydot, [], 1);
        end
        
        function jacobian = get_jacobian_sparsity(obj)
            % gets the jacobian sparsity matrix for this model
            % (this speeds up computation)
            
            dNextCell = obj.get_dNextCell();
            dPrevCell = obj.get_dPrevCell();
            dCurrCell = obj.get_dCurrCell();
            
            S = zeros(obj.numSpecies*obj.N, obj.numSpecies*obj.N);            

            % cell 1
            S(1 : obj.numSpecies , 1:obj.numSpecies               ) = dCurrCell;
            S(1 : obj.numSpecies , obj.numSpecies+1 : 2*obj.numSpecies) = dNextCell;
            
            % cells 2 through N-1
            for i = obj.numSpecies+1 : obj.numSpecies : obj.numSpecies*(obj.N-1)
                S(i : i+(obj.numSpecies-1)  , i - obj.numSpecies    : i - 1                   )  = dPrevCell;
                S(i : i+(obj.numSpecies-1)  , i                     : i + (obj.numSpecies-1)  )  = dCurrCell;
                S(i : i+(obj.numSpecies-1)  , i + obj.numSpecies    : i + (2*obj.numSpecies-1)    )  = dNextCell;
            end
            
            % cell N
            S(obj.numSpecies*(obj.N-1)+1 : obj.numSpecies*obj.N     , obj.numSpecies*(obj.N-2)+1    : obj.numSpecies*(obj.N-1)) = dPrevCell;
            S(obj.numSpecies*(obj.N-1)+1 : obj.numSpecies*obj.N     , obj.numSpecies*(obj.N-1)+1    : obj.numSpecies*obj.N  ) = dCurrCell;
            
            % output is sparsified S.  
            jacobian = sparse(S);  
        end
        
        function species = get_Ss_species_1D(obj)
            % returns indices of species that should be at steady-state in
            % the 1D model
            species = (1 : obj.numSpecies * obj.N);
        end
        
        function species = get_Ss_species_0D(obj)
            % returns indices of species that should be at steady-state in
            % the 0D model 
            species = (1 : obj.numSpecies);
        end
        
        function initial_conditions = get_initial_conditions_0D(obj)
            % returns initial conditions for the 0D preequilibration step
            % assume that PTCH is in the steady state, and none has been
            % converted to C yet
            initial_conditions = zeros(1, obj.numSpecies - 1);
            % ^ replace in subclasses
        end
        
        function value = get_dPrevCell(obj)
            value = zeros(obj.numSpecies);
            value(1, 1) = 1;  % set dh(i)/dh(i-1) equal to 1
        end
        
        function value = get_dNextCell(obj)
            value = zeros(obj.numSpecies);
            value(1, 1) = 1;  % set dh(i)/dh(i-1) equal to 1
        end
        
        function value = get_dCurrCell(obj)
            % dependence on same cell (dy(i) / dy(i))
            value = zeros(obj.numSpecies);
        end
    end
end