% ----------------------------------------------------------------------- %
% Function MDOMOPSO performs a Multi-Objective Particle Swarm Optimization%
% over continous functions including the storing of Multiple Design       %
% Options                                                                 %
%                                                                         %
%   Input parameters:                                                     %
%       - fun: Anonymous multi-obj function to minimize.                  %
%           The function expected type is [objs, varoutput] = fun(x)      %
%               -objs: is a vector of the objective functions to minimize %
%               -varoutput [opt]: vector of auxiliary values to be stored %
%       - nVar:    Number of variables.                                   %
%       - var_min: Vector that indicates the minimum values of the search %
%                                           space in each dimension.      %
%       - var_max: Same than 'var_min' with the maxima.                   %
%       - options: Struct containing optional values for the algorithm    %
%           Admitted parameters for the optimization algorithm            %
%           * Np:        Number of particles.                             %
%           * Nr:        Repository size (in particles).                  %
%           * maxgen:    Maximum number of generations.                   %
%           * W:         Inertia coefficient.                             %
%           * C1:        Personal confidence factor.                      %
%           * C2:        Swarm confidence factor.                         %
%           * ngrid:     Number of hypercubes in each dimension.          %
%           * maxvel:    Maximum velocity (search space percentage)       %
%           * u_mut:     Uniform mutation percentage.                     %
%           * aux_output: true to store auxiliary output.                 %
%           * vectorized: true if the function is vectorizable.           %
%           * parallel:  if true the algorithm exploits parallelization.  %
%                       (applies when vectorized is false)                %
%           * stall_iter: (NOT USED) Stall iterations after which the     %
%                           methodology converges                         %
%           * additional_iter: additional iterations to perform after     %
%                        stopping criterion has been reached              %
%           * tolParetoChange: tolerance on the Pareto front change       %
%           * verbose: if true, iteration status is written               %
%           * MDO_calc: if true, the procedure calculates the MDO points  %
%           * MDO_tol: tolerance of the MDO points                        %
%                     if 2 objective functions are considered, then,      %
%                     all MDO points within MDO_tol with respect to the   %
%                     translated convex hull of Pareto frontier is        %
%                     selected, otherwise, all MDOs within a MDO_tol with %
%                     respect to the final Pareto frontier are selected   %
%                                                                         %
%   Output:                                                               %
%       REP structure cointaining the Pareto frontier and the history     %
%       * REP.history: structure of the history of the evaluated points   %
%       * REP.history.pos: vector of the evaluated particles' position    % 
%       * REP.history.fit: vector of the fitness value for each particle  % 
%       * REP.history.aux: auxiliary quantities for each particleition    %
%       * REP.history.Spread: Spread by iteration                         %
%       * REP.history.avgDistance: avgDistance                            %
%       * REP.pos: position of the particles of the final Pareto frontier %
%       * REP.fit: fitness value for of the particles in REP.pos          %
%                                                                         %
%       MDO structure cointaining the MDO points                          %
%       * MDO.pos: vector of the MDO position                             % 
%       * MDO.fit: vector of the MDOfitness value                         % 
%       * MDO.aux: auxiliary quantities for each MDO                      %
%       * MDO.selected_ids: ids of the MDOs                            %
%                                                                         %
% ----------------------------------------------------------------------- %
%   For an example of use, run 'example.m'.                               %
% ----------------------------------------------------------------------- %
%   Author:  Davide Fioriti                                               %
%                                                                         %
%   Code modified from code by Victor Martinez Cagigal Copyright (c) 2017 %
%   Original code available at:                                           %
%   https://it.mathworks.com/matlabcentral/fileexchange/                  %
%        62074-multi-objective-particle-swarm-optimization-mopso          %
%                                                                         %
%   Date:    2/11/2020                                                    %
%   E-mail:  fioritidavidesubs (at) gmail (dot) com                       %
%   Version: 1                                                            %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%    [1]Fioriti, D., Lutzemberger, G., Poli, D., Duenas-Martinez, P.,     %
%       Micangeli, A., Coupling economic multi-objective optimization     %
%       and multiple design options: a business-oriented approach to      %
%       optimize an off-grid hybrid microgrid, International Journal of   %
%       Electrical Power and Energy Systems, 2021                         %
%                                                                         %
%    [2]Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
%                                                                         %
%    [3]Sierra, M. R., & Coello, C. A. C. (2005, March). Improving PSO-   %
%       based multi-objective optimization using crowding, mutation and ?-%
%       dominance. In International Conference on Evolutionary Multi-Crite%
%       rion Optimization (pp. 505-519). Springer Berlin Heidelberg.      %
% ----------------------------------------------------------------------- %
function [REP, MDO] = MDOMOPSO(fun, nVar, var_min, var_max, options)

    % Parameters
    params.Np = 200;        % Population size
    params.Nr = 100;        % Repository size
    params.maxgen = 100;    % Maximum number of generations
    params.W = 0.4;         % Inertia weight
    params.C1 = 2;          % Individual confidence factor
    params.C2 = 2;          % Swarm confidence factor
    params.ngrid = 20;      % Number of grids in each dimension
    params.maxvel = 5;      % Maxmium vel in percentage
    params.u_mut = 0.5;     % Uniform mutation percentage
    params.vectorized = false; % True if the function is vectorized
    params.aux_output = true; % True if additional output should be stored
    params.parallel = false; % When true, parallel computing is used
    params.stall_iter = 10;   % NOT USED Stall iterations on which to calculate the spread
    params.additional_iter = 0;   % Additional iterations to compute after a stopping criterion has been reached
    params.tolParetoChange = 1e-5; % Tolerance on the Pareto change that if kept for stall_iter iterations, the algorithm stops
    params.verbose = true; % If true, the iteration status is written
    params.MDO_calc = true; % If true, it calculates the points within MDO_tolerance with respect to the pareto frontier
    params.MDO_tol = 0.05; % Tolerance of the MDO points with respect to the pareto frontier
    
    if nargin < 4
        error("missing input variables")
    elseif nargin >4
        S=fieldnames(options);
        for i = 1:numel(S)
            if isfield(params, S{i})
                params.(S{i}) = options.(S{i});
            else
                error("Option not recognized")
            end
        end
    end
    
    % Initialization
    var_min = var_min(:);
    var_max = var_max(:);
    
    POS = repmat((var_max-var_min)',params.Np,1).*rand(params.Np,nVar) + repmat(var_min',params.Np,1);
    VEL = zeros(params.Np,nVar);
    
    % Evaluate the population
    [POS_fit, aux] = fun_vect(fun, POS, params.Np, [], [], params.aux_output, params.vectorized, params.parallel);

    % Evaluate the number of outputs of the function
    n_objs = size(POS_fit, 2);
    
    % Evaluate the number of auxiliary data
    n_aux = size(aux, 2);
    
    % Initialization of history
    REP.history.pos = zeros(params.maxgen*params.Np, nVar); % history vector population position
    REP.history.fit = zeros(params.maxgen*params.Np, n_objs); % history vector fitness function
    REP.history.aux = zeros(params.maxgen*params.Np, n_aux); % history vector fitness function
        
    if size(POS,1) ~= size(POS_fit,1)
        warning(['The objective function is badly programmed. It is not returning' ...
            'a value for each particle, please check it.']);
    end
    PBEST    = POS;
    PBEST_fit= POS_fit;
    DOMINATED= checkDomination(POS_fit);
    REP.pos  = POS(~DOMINATED,:);
    REP.fit = POS_fit(~DOMINATED,:);
    REP      = updateGrid(REP,params.ngrid);
    maxvel   = (var_max-var_min).*params.maxvel./100;
    gen      = 1;
    
    REP.history.pos(1:params.Np, :) = POS;
    REP.history.fit(1:params.Np, :) = POS_fit;
    REP.history.aux(1:params.Np, :) = aux;
    Spread = zeros(params.maxgen, 1); %Spread measure for the convergence criteria
    avgDistance = zeros(params.maxgen, 1); %Normalized distance measure of the pareto front
    
    % Plotting and verbose
    if(size(POS_fit,2)==2) && params.verbose
        h_fig = figure(1);
        h_par = plot(POS_fit(:,1),POS_fit(:,2),'or'); hold on;
        h_rep = plot(REP.fit(:,1),REP.fit(:,2),'ok'); hold on;
        try
            set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
            axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                  min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            grid on; xlabel('f1'); ylabel('f2');
        end
        drawnow;
    end
    if(size(POS_fit,2)==3) && params.verbose
        h_fig = figure(1);
        h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
        h_rep = plot3(REP.fit(:,1),REP.fit(:,2),REP.fit(:,3),'ok'); hold on;
        try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
        end
        grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
        drawnow;
        axis square;
    end
    if params.verbose
        fprintf('\n\n%10s %10s %8s %8s\n', 'Iteration', 'Pareto p.', 'Spread', 'AvgDistance');
        fprintf('%10d %10d %8.4f %8.4f\n', gen, size(REP.pos,1), Spread(gen), avgDistance(gen));
    end
    
    % Main MPSO loop
    stopCondition = false;
    stop_criterion_reached = false;
    convergence_iter = Inf;
    while ~stopCondition 
        gen = gen + 1;
        
        % Select leader
        h = selectLeader(REP);
        
        % Update speeds and positions
        VEL = params.W.*VEL + params.C1*rand(params.Np,nVar).*(PBEST-POS) ...
                     + params.C2*rand(params.Np,nVar).*(repmat(REP.pos(h,:),params.Np,1)-POS);
        POS = POS + VEL;
        
        % Perform mutation
        POS = mutation(POS,gen,params.maxgen,params.Np,var_max,var_min,nVar,params.u_mut);
        
        % Check boundaries
        [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min);       
        
        % Evaluate the population
        [POS_fit, aux] = fun_vect(fun, POS, params.Np, n_objs, n_aux, params.aux_output, params.vectorized, params.parallel);
        
        %update history
        REP.history.pos((gen-1)*params.Np + (1:params.Np), :) = POS; % update history of population position
        REP.history.fit((gen-1)*params.Np + (1:params.Np), :) = POS_fit; % update history of fitness function
        REP.history.aux((gen-1)*params.Np + (1:params.Np), :) = aux; % update history of auxiliary data
        
        %store the previous repository of the fitness functions
        pre_pos_fit = REP.fit;
        
        % Update the repository
        REP = updateRepository(REP,POS,POS_fit,params.ngrid);
        if(size(REP.pos,1)>params.Nr)
            crowding_distances = calculate_crowding_distances(REP.fit);
            REP = deleteFromRepository(REP, crowding_distances, size(REP.pos,1)-params.Nr,params.ngrid);
        end
        crowding_distances = calculate_crowding_distances(REP.fit);
        Spread(gen) = calculate_spread(crowding_distances, REP.fit, pre_pos_fit);
        avgDistance(gen) = calculate_avgDistance(crowding_distances);
        
        % Update the best positions found so far for each particle
        pos_best = dominates(POS_fit, PBEST_fit);
        best_pos = ~dominates(PBEST_fit, POS_fit);
        best_pos(rand(params.Np,1)>=0.5) = 0;
        if(sum(pos_best)>1)
            PBEST_fit(pos_best,:) = POS_fit(pos_best,:);
            PBEST(pos_best,:) = POS(pos_best,:);
        end
        if(sum(best_pos)>1)
            PBEST_fit(best_pos,:) = POS_fit(best_pos,:);
            PBEST(best_pos,:) = POS(best_pos,:);
        end
        
        % Plotting and verbose
        if(size(POS_fit,2)==2) && params.verbose
            figure(h_fig); delete(h_par); delete(h_rep);
            h_par = plot(POS_fit(:,1),POS_fit(:,2),'or'); hold on;
            h_rep = plot(REP.fit(:,1),REP.fit(:,2),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2))]);
            end
            grid on; xlabel('f1'); ylabel('f2');
            drawnow;
            axis square;
        end
        if(size(POS_fit,2)==3) && params.verbose
            figure(h_fig); delete(h_par); delete(h_rep); 
            h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
            h_rep = plot3(REP.fit(:,1),REP.fit(:,2),REP.fit(:,3),'ok'); hold on;
            try
                set(gca,'xtick',REP.hypercube_limits(:,1)','ytick',REP.hypercube_limits(:,2)','ztick',REP.hypercube_limits(:,3)');
                axis([min(REP.hypercube_limits(:,1)) max(REP.hypercube_limits(:,1)) ...
                      min(REP.hypercube_limits(:,2)) max(REP.hypercube_limits(:,2)) ...
                      min(REP.hypercube_limits(:,3)) max(REP.hypercube_limits(:,3))]);
            end
            grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
            drawnow;
            axis square;
        end
        if params.verbose
            fprintf('%10d %10d %8.4f %8.4f\n', gen, size(REP.pos,1), Spread(gen), avgDistance(gen));
        end
        
        % Update generation and check for termination
        if(gen>=params.maxgen)
            stopCondition = true;
        elseif stop_criterion_reached
            if (gen - convergence_iter >= params.additional_iter)
                stopCondition = true;
            end
        elseif gen > params.stall_iter % stall iterations not implemented yet
            if abs(Spread(gen-1) - Spread(gen)) <= params.tolParetoChange*max(1,Spread(gen-1)) || ...
                    abs(avgDistance(gen-1) - avgDistance(gen)) <= params.tolParetoChange*max(1,avgDistance(gen-1))
                stop_criterion_reached = true;
                convergence_iter = gen;
                if params.additional_iter == 0
                    stopCondition = true;
                end
                if params.verbose
                    fprintf('%10d Criteria reached\n', gen);
                end
            end
        end
    end
    
    % Remove unused data in the history
    REP.history.pos(gen*params.Np+1:end, :) = [];
    REP.history.fit(gen*params.Np+1:end, :) = [];
    REP.history.aux(gen*params.Np+1:end, :) = [];
    REP.history.Spread = Spread(1:gen, :);
    REP.history.avgDistance = avgDistance(1:gen, :);
    hold off;
    
    %initialization of MDO structure
    MDO = struct();
    if params.MDO_calc
        % Calculate only if MDO_calc is enabled
        
        if n_objs == 2
            % 2 objective functions
            
            [~, ids] = sort(REP.fit(:,1));
            fval_PARETO_sort = REP.fit(ids, :);
            
            if numel(params.MDO_tol) == 1
                params.MDO_tol = ones(1, n_objs) * params.MDO_tol;
            else
                params.MDO_tol = params.MDO_tol(:)';
            end
            
            max_1 = max(fval_PARETO_sort(:,1))*(1+params.MDO_tol(1));
            max_2 = max(fval_PARETO_sort(:,2))*(1+params.MDO_tol(2));

            polygon = [0,0; 0, max_2; fval_PARETO_sort(:, :).*(1+params.MDO_tol); max_1,0; 0,0];

            MDO.selected_ids = inpolygon(REP.history.fit(:,1), REP.history.fit(:,2), polygon(:,1), polygon(:,2));
        else
            % multiple objective functions
            
            MDO_ids_raw = arrayfun(@(id) find(all(REP.history.fit < (REP.fit(id, :) + abs(REP.fit(id, :)) .* params.MDO_tol), 2)), 1:size(REP.pos, 1), 'UniformOutput', false);
            MDO_ids = [];
            for id = 1:size(REP.pos, 1)
                MDO_ids = [MDO_ids; MDO_ids_raw{id}];
            end
            
            MDO.selected_ids = sort(unique(MDO_ids));
        end
        
        MDO.fit = REP.history.fit(MDO.selected_ids,:);
        MDO.pos = REP.history.pos(MDO.selected_ids,:);
        MDO.aux = REP.history.aux(MDO.selected_ids,:);
    end
end

% This is an itnermediate function to allow vectorization of non-vectorized
%       function and to manage the auxiliary output functions
function [POS_fit, auxiliary] = fun_vect(fun, x, Np, n_objs, n_aux, aux_output, vectorized, parallel_for)
    if vectorized
        if aux_output
            [POS_fit, auxiliary] = fun(x);
        else
            [POS_fit] = fun(x);
            auxiliary = zeros(Np, 0);
        end
    else
        if isempty(n_objs) || n_objs < 0 || isempty(n_aux) || n_aux < 0
            if aux_output
                [POS_temp, aux_temp] = fun(x(1, :));
            else
                POS_temp = fun(x(1, :));
                aux_temp = [];
            end
            n_objs = numel(POS_temp);
            n_aux = numel(aux_temp);
        end
        POS_fit = zeros(Np, n_objs);
        auxiliary = zeros(Np, n_aux);
        if aux_output
            if parallel_for
                parfor i=1:Np
                    [POS_fit(i, :), auxiliary(i, :)] = fun(x(i, :));
                end
            else
                for i=1:Np
                    [POS_fit(i, :), auxiliary(i, :)] = fun(x(i, :));
                end
            end
        else
            if parallel_for
                parfor i=1:Np
                    POS_fit(i, :) = fun(x(i, :));
                end
            else
                for i=1:Np
                    POS_fit(i, :) = fun(x(i, :));
                end
            end
        end
    end
end

% Function that updates the repository given a new population and its
% fitness
function REP = updateRepository(REP,POS,POS_fit,ngrid)
    % Domination between particles
    DOMINATED  = checkDomination(POS_fit);
    REP.pos    = [REP.pos; POS(~DOMINATED,:)];
    REP.fit= [REP.fit; POS_fit(~DOMINATED,:)];
    % Domination between nondominated particles and the last repository
    DOMINATED  = checkDomination(REP.fit);
    REP.fit= REP.fit(~DOMINATED,:);
    REP.pos    = REP.pos(~DOMINATED,:);
    % Updating the grid
    REP        = updateGrid(REP,ngrid);
end

% Function that corrects the positions and velocities of the particles that
% exceed the boundaries
function [POS,VEL] = checkBoundaries(POS,VEL,maxvel,var_max,var_min)
    % Useful matrices
    Np = size(POS,1);
    MAXLIM   = repmat(var_max(:)',Np,1);
    MINLIM   = repmat(var_min(:)',Np,1);
    MAXVEL   = repmat(maxvel(:)',Np,1);
    MINVEL   = repmat(-maxvel(:)',Np,1);
    
    % Correct positions and velocities
    VEL(VEL>MAXVEL) = MAXVEL(VEL>MAXVEL);
    VEL(VEL<MINVEL) = MINVEL(VEL<MINVEL);
    VEL(POS>MAXLIM) = (-1).*VEL(POS>MAXLIM);
    POS(POS>MAXLIM) = MAXLIM(POS>MAXLIM);
    VEL(POS<MINLIM) = (-1).*VEL(POS<MINLIM);
    POS(POS<MINLIM) = MINLIM(POS<MINLIM);
end

% Function for checking the domination between the population. It
% returns a vector that indicates if each particle is dominated (1) or not
function dom_vector = checkDomination(fitness)
    Np = size(fitness,1);
    dom_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
    d = all(x<=y,2) & any(x<y,2);
end

% Function that updates the hypercube grid, the hypercube where belongs
% each particle and its quality based on the number of particles inside it
function REP = updateGrid(REP,ngrid)
    % Computing the limits of each hypercube
    ndim = size(REP.fit,2);
    REP.hypercube_limits = zeros(ngrid+1,ndim);
    for dim = 1:1:ndim
        REP.hypercube_limits(:,dim) = linspace(min(REP.fit(:,dim)),max(REP.fit(:,dim)),ngrid+1)';
    end
    
    % Computing where belongs each particle
    npar = size(REP.fit,1);
    REP.grid_idx = zeros(npar,1);
    REP.grid_subidx = zeros(npar,ndim);
    for n = 1:1:npar
        idnames = [];
        for d = 1:1:ndim
            REP.grid_subidx(n,d) = find(REP.fit(n,d)<=REP.hypercube_limits(:,d)',1,'first')-1;
            if(REP.grid_subidx(n,d)==0), REP.grid_subidx(n,d) = 1; end
            idnames = [idnames ',' num2str(REP.grid_subidx(n,d))];
        end
        REP.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames ');']);
    end
    
    % Quality based on the number of particles in each hypercube
    REP.quality = zeros(ngrid,2);
    ids = unique(REP.grid_idx);
    for i = 1:length(ids)
        REP.quality(i,1) = ids(i);  % First, the hypercube's identifier
        REP.quality(i,2) = 10/sum(REP.grid_idx==ids(i)); % Next, its quality
    end
end

% Function that selects the leader performing a roulette wheel selection
% based on the quality of each hypercube
function selected = selectLeader(REP)
    % Roulette wheel
    prob    = cumsum(REP.quality(:,2));     % Cumulated probs
    sel_hyp = REP.quality(find(rand(1,1)*max(prob)<=prob,1,'first'),1); % Selected hypercube
    
    % Select the index leader as a random selection inside that hypercube
    idx      = 1:1:length(REP.grid_idx);
    selected = idx(REP.grid_idx==sel_hyp);
    selected = selected(randi(length(selected)));
end

% Calculate crowding distances
function crowding_distances = calculate_crowding_distances(data)
    % Compute the crowding distances
    crowding_distances = zeros(size(data,1),1);
    for m = 1:1:size(data,2)
        [m_fit,idx] = sort(data(:,m),'ascend');
        m_up     = [m_fit(2:end); Inf];
        m_down   = [Inf; m_fit(1:end-1)];
        distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
        [~,idx]  = sort(idx,'ascend');
        crowding_distances = crowding_distances + distance(idx);
    end
    crowding_distances(isnan(crowding_distances)) = Inf;
end

% Calculate Spread
function Spread = calculate_spread(crowding_distances, score, score_old)
    is_finite = isfinite(crowding_distances);
    if isempty(is_finite)
        Spread = 0;
    else
        mean_crowd = mean(crowding_distances(is_finite));
        std_crowd = std(crowding_distances(is_finite));
        Q = size(score, 2);

        mu = 0;
        for k = 1:size(score, 2)
            [min_new, min_pos_new] = min(score(:, k));
            min_new_row = score(min_pos_new, :);
            [min_old, min_pos_old] = min(score_old(:, k));
            min_old_row = score_old(min_pos_old, :);

            mu = mu + norm(min_new_row - min_old_row);
        end
        
        if not(isfinite(mu))
            mu = 0;
        end
        
        if mu == 0 && mean_crowd == 0
            Spread = 0;
        else
            Spread = (mu + std_crowd)/(mu + Q * mean_crowd);
        end
    end
end

%Calculate average distance
function avgDistance = calculate_avgDistance(crowding_distances)
    crowding_distances = crowding_distances(isfinite(crowding_distances));
    if isempty(crowding_distances)
        avgDistance = 0;
    else
        %avgDistance = norm(crowding_distances - mean(crowding_distances))/sqrt(length(crowding_distances));
        avgDistance = sqrt(sum(crowding_distances.^2)/length(crowding_distances));
    end
end

% Function that deletes an excess of particles inside the repository using
% crowding distances
function REP = deleteFromRepository(REP,crowding_distances, n_extra,ngrid)    
    % Delete the extra particles with the smallest crowding distances
    [~,del_idx] = sort(crowding_distances,'ascend');
    del_idx = del_idx(1:n_extra);
    REP.pos(del_idx,:) = [];
    REP.fit(del_idx,:) = [];
    REP = updateGrid(REP,ngrid); 
end

% Function that performs the mutation of the particles depending on the
% current generation
function POS = mutation(POS,gen,maxgen,Np,var_max,var_min,nVar,u_mut)
    % Sub-divide the swarm in three parts [2]
    fract     = Np/3 - floor(Np/3);
    if(fract<0.5), sub_sizes =[ceil(Np/3) round(Np/3) round(Np/3)];
    else           sub_sizes =[round(Np/3) round(Np/3) floor(Np/3)];
    end
    cum_sizes = cumsum(sub_sizes);
    
    % First part: no mutation
    % Second part: uniform mutation
    nmut = round(u_mut*sub_sizes(2));
    if(nmut>0)
        idx = cum_sizes(1) + randperm(sub_sizes(2),nmut);
        POS(idx,:) = repmat((var_max-var_min)',nmut,1).*rand(nmut,nVar) + repmat(var_min',nmut,1);
    end
    
    % Third part: non-uniform mutation
    per_mut = (1-gen/maxgen)^(5*nVar);     % Percentage of mutation
    nmut    = round(per_mut*sub_sizes(3));
    if(nmut>0)
        idx = cum_sizes(2) + randperm(sub_sizes(3),nmut);
        POS(idx,:) = repmat((var_max-var_min)',nmut,1).*rand(nmut,nVar) + repmat(var_min',nmut,1);
    end
end