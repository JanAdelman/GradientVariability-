function readout_precision

% options
simulate = true; % if false, plot results from saved files instead of generating new data
write = true; % when creating new data, should the results be written to output files?
fitcosh = true; % fit an exp or cosh to the gradients
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 100; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 4.9; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 19.26; % mean gradient length [µm]
mu_a = (diameter/2)^2*pi; %mean cell area 
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain
%CV = [0.01:0.01:0.09 0.1:0.05:0.95 1:0.5:10]'; % coefficient of variation of the kinetic parameters
CV = [0.3];
%plot_CV = CV(14); % plot the gradients for this CV value
CV_area = 0.5;
% analytical deterministic solution
nc = ncS + ncP; % total number of cells
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

%CVfun = @(x) nanstd(x) / nanmean(x);
CVfun = @(x) nanstd(x) ./ nanmean(x);
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));
Stdfun = @(x) nanstd(x);

fitopt = statset('TolFun', tol, 'TolX', tol);

close all
%f1 = figure('Name', 'Individual gradients', 'Position', [0 0 2000 800]);
%names = {'p', 'd', 'D', 'all'};
names = {'p'};

% get the readout positions for the gradient 
readout_pos = [0:LP/100:LP];
%% vary molecular noise

f2 = figure('Name', 'Spread of Variability in Concentration for different Positions', 'Position', [0 0 2000 800]);
%f3 = figure('Name', 'Dependency of gradient variability on molecular noise', 'Position', [0 0 2000 800]);
                    
% k = p, d, D, and all three together
for k = 1:numel(names)
    filename = ['read_out_precision_average_conc_' names{k} '.csv'];
    if names{k} == 'D'
        filename = 'read_out_precision_average_conc_diff.csv';
    end
        
    if simulate
        
        conc_average = NaN(length(CV), length(readout_pos));
        conc_SE_average = NaN(length(CV), length(readout_pos));
        conc_CV_average = NaN(length(CV), length(readout_pos));
        conc_CV_SE_average = NaN(length(CV), length(readout_pos));
        conc_random = NaN(length(CV), length(readout_pos));
        conc_SE_random = NaN(length(CV), length(readout_pos));
        conc_CV_SE_random = NaN(length(CV), length(readout_pos));
        conc_CV_random = NaN(length(CV), length(readout_pos));
        conc_cilium = NaN(length(CV), length(readout_pos));
        conc_SE_cilium = NaN(length(CV), length(readout_pos));
        conc_CV_cilium = NaN(length(CV), length(readout_pos));
        conc_CV_SE_cilium = NaN(length(CV), length(readout_pos));
        conc_std_average = NaN(length(CV), length(readout_pos));
        conc_std_random = NaN(length(CV), length(readout_pos));
        conc_std_cilium = NaN(length(CV), length(readout_pos));
               
        % loop over variabilities
        for i = 1:length(CV)
            
            % loop over several independent runs
            conc_per_iteration_average = NaN(nruns, length(readout_pos));
            conc_per_iteration_random = NaN(nruns, length(readout_pos));
            conc_per_iteration_cilium = NaN(nruns, length(readout_pos));
            
            for j = 1:nruns
                
                 % initialise arrays for variable cell size          
                l_s_temp = [];
                l_s = [];
                            
                l_p_temp = [];
                l_p = [];
                
                % ======================================= %
                % Area computations for the Source Domain %
                % ======================================= %

                while sum(l_s_temp) < LS 
                    
                    rand_area = random(logndist(mu_a, mu_a * CV_area), 1, 1);
                    length_s_sampled = 2*sqrt((rand_area)/pi);
                    
                    
                    if rand_area >= LS
                        
                        l_s_temp = [l_s_temp, length_s_sampled];
                        
                        
                    end
                    
                    l_s_temp = [l_s_temp, length_s_sampled];
                    
                end
                          
                sum_s_upper = sum(l_s_temp);
                sum_s_lower = sum(l_s_temp(1:end-1));
                
                % edge case, where the source consists only of one cell 
                if length(l_s_temp) == 1
                
                    l_s = l_s_temp;
                    
                else    
                    if abs(sum_s_upper - LS) < abs(sum_s_lower - LS)
                        l_s = l_s_temp;
                    else
                        l_s = l_s_temp(1:end-1);
                    end 
                end
                
                              
                % from the area calculate the diameter of each cell 
                d_s_normalised = l_s / sum(l_s);
                l_s = d_s_normalised * LS;
                l_s = cumsum(l_s);                          
                l_s = fliplr(l_s);
                
                % =========================================== %
                % Area computations for the Patterning Domain %
                % =========================================== %
                
                 % add cells as long as the overall arrea is not surpassed
                 while sum(l_p_temp) < LP 
                    
                    rand_area = random(logndist(mu_a, mu_a * CV_area), 1, 1);
                    length_p_sampled = 2*sqrt((rand_area)/pi);
                    
                    if rand_area >= LP
                        
                        l_p_temp = [l_p_temp, length_p_sampled];
                    end
                    
                    
                    l_p_temp = [l_p_temp, length_p_sampled];
                    
                 end
                
                % add one more cell and check if the area is closer to the
                % desired area compaed to the case wihout this last cell
                sum_p_upper = sum(l_p_temp);
                sum_p_lower = sum(l_p_temp(1:end-1));
                        
                if abs(sum_p_upper - LP) < abs(sum_p_lower - LP)
                    l_p = l_p_temp;
                else
                    l_p = l_p_temp(1:end-1);
                end              
        
                % calculate the normalied diameter for each cell               
                d_p_normalised = l_p /sum(l_p);
                l_p = d_p_normalised * LP;
                l_p = cumsum(l_p);
                
                 % initialise the solver
                x0 = [];
                x0 = [x0, -l_s, 0, l_p];
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
                
                nc = length(l_p) + length(l_s);
                ncS = length(l_s);
                ncP = length(l_p);
                
                options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);
                
                 % default: all parameters constant
                p = mu_p * ones(nc, 1);
                d = mu_d * ones(nc, 1);
                D = mu_D * ones(nc, 1);

                % draw random kinetic parameters for each cell
                if k == 1 || k == 4
                    p = random(logndist(mu_p, mu_p * CV(i)), nc, 1);
                end
                if k == 2 || k == 4
                    d = random(logndist(mu_d, mu_d * CV(i)), nc, 1);
                end
                if k == 3 || k == 4
                    D = random(logndist(mu_D, mu_D * CV(i)), nc, 1);
                end
                p = random(logndist(mu_p, mu_p * CV(i)), nc, 1);
                d = random(logndist(mu_d, mu_d * CV(i)), nc, 1); 
                D = random(logndist(mu_D, mu_D * CV(i)), nc, 1);
                
                % get initial solution 
                sol0 = bvpinit(x0, @y0);
                
                % solve the equation
                sol = bvp4c(@odefun, @bcfun, sol0, options);
                
                % initialise the start location at the beginning of ther
                % patterning domain 
                cell_beginning = 0;
                
                % array to store numerical integration results
                y_sol_average = [];
                
                % array to store random conc. from a cell
                y_sol_random = [];
                
                % array to store conc. closeste to mid point
                y_sol_cilium = [];
                
                % use the beginning of each cell as the x-coordinate 
                x_sol = [0, l_p(1:end-1)];
                               
                % calculate the middle point of each cell 
                middle = (l_p - x_sol) /2;
                mid_point = x_sol + middle;
                              
                % loop through the  cell and extract the solutions for each
                % cell separately. Then use the trapezoid method to
                % numerically integrate and find the mean morphogen
                % gradient concentration over one cell.
                
                for cell_loc = 1:ncP
                    
                    % set the upper interval as the end of a cell 
                    cell_end = l_p(cell_loc);
            
                    % define interval where to extract solutions 
                    logical_indexes = (sol.x <= cell_end) & (sol.x >= cell_beginning);
                    % extract indices of the desired solutions 
                    interval = find(logical_indexes);
                    
                    % get lenght of the cell for normalisation
                    cell_length = cell_end - cell_beginning;
                                                      
                    % set the lower interval for the next iteration as the
                    % current end of the cell 
                    cell_beginning = cell_end;
                    
                    % get the x and y solution 
                    X = sol.x(interval);
                    Y = sol.y(1, interval);
                    
                    % get the x value closest to the mean cell area
                    
                    [min_distance, index] = min( abs(X - mid_point(cell_loc)));
                  
                    cilium_sol = Y(index);
                    
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_cilium = [y_sol_cilium, cilium_sol];
                    
                    % get a concentration randomly from a cell 
                    rand_conc = randsample(Y,1);
                    
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_random = [y_sol_random, rand_conc];
                    
                                     
                    % get the average concentration per cell (thus
                    % normalised by cell length) 
                    trapz_sol = trapz(X,Y)/cell_length;
  
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_average = [y_sol_average, trapz_sol];
                    
                                      
                end 
                
                piecewise_const_average = [];
                piecewise_const_random = [];
                piecewise_const_cilium = [];
                
                % loop through readout position to get the concenctration
                % c(x) at each x value. (Piecewise constant function)
                for idx = 1:length(readout_pos)
                    
                    % if at tthe end, set the index to the position of the
                    % last cell 
                    if readout_pos(idx) >= x_sol(end)
                        B = length(x_sol);
                        piecewise_const_average = [piecewise_const_average, y_sol_average(B)];
                        piecewise_const_random = [piecewise_const_random, y_sol_random(B)];
                        piecewise_const_cilium = [piecewise_const_cilium, y_sol_cilium(B)];
                        
                    end 
                    
                    % find indices where the readout position is smaller
                    % than the x value defined by a cell. The Index one to
                    % the left defines the value of of the concentration
                    % c(x). 
                    
                    B = find(readout_pos(idx) < x_sol);
                    B = min(B)-1;
                    
                    piecewise_const_average = [piecewise_const_average, y_sol_average(B)]; 
                    piecewise_const_random = [piecewise_const_random, y_sol_random(B)];
                    piecewise_const_cilium = [piecewise_const_cilium, y_sol_cilium(B)];
                                        
                end
                % get all concenctrations for one iteration at each
                % location x. 
                conc_per_iteration_average(j, :) = piecewise_const_average;
                conc_per_iteration_random(j, :) = piecewise_const_random;
                conc_per_iteration_cilium(j, :) = piecewise_const_cilium;
                
            end
            
            % Caculate summary statistics for the nruns
            % calculate the mean concentration at each position        
            conc_average(i, :) = mean(conc_per_iteration_average);
            conc_random(i, :) = mean(conc_per_iteration_random);
            conc_cilium(i, :) = mean(conc_per_iteration_cilium);
                                
            % calculate the standard error at each position               
            conc_SE_average(i, :) = SEfun(conc_per_iteration_average);
            conc_SE_random(i, :) = SEfun(conc_per_iteration_random);
            conc_SE_cilium(i, :) = SEfun(conc_per_iteration_cilium);
            
            % calculate standard deviation at each position
            conc_std_average(i, :) = Stdfun(conc_per_iteration_average);
            conc_std_random(i, :) = Stdfun(conc_per_iteration_random);
            conc_std_cilium(i, :) = Stdfun(conc_per_iteration_cilium);
            
            % calculate the coefficient of variation at each position 
            conc_CV_average(i, :) = CVfun((conc_per_iteration_average));
            conc_CV_random(i, :) = CVfun((conc_per_iteration_random));
            conc_CV_cilium(i, :) = CVfun((conc_per_iteration_cilium));
            
            % calculate the SE for the coefficient of variation 
            conc_CV_SE_average(i, :) =  std(bootstrp(nboot, CVfun, conc_per_iteration_average));
            conc_CV_SE_random(i, :) =  std(bootstrp(nboot, CVfun, conc_per_iteration_random));
            conc_CV_SE_cilium(i, :) =  std(bootstrp(nboot, CVfun, conc_per_iteration_cilium));
            
        end
              
        % write data
        if write
            writetable(table(conc_average, conc_SE_average, conc_CV_average, conc_CV_SE_average), filename);
        end
    else
        % read data
      
    end
   
    % plot the relationship between CV_k and lambda
    figure(f2)
    
    subplot(3,1,1)
    errorbar(readout_pos, conc_cilium(i, :), conc_SE_cilium(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'g')   
    hold on
    errorbar(readout_pos, conc_random(i, :), conc_SE_random(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'r')
    errorbar(readout_pos, conc_average(i, :), conc_SE_average(i, :), 'bo', 'LineWidth', LineWidth)   
    xlabel(['Position x [µm]'])
    ylabel('µ_{C(x)}')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'linear')
    grid on
    
    subplot(3,1,2);
    %errorbar(readout_pos, conc_std_cilium(i, :), conc_SE_cilium(i, :),'bo', 'LineWidth', LineWidth, 'Color', 'g')   
    plot(readout_pos, conc_std_cilium(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'g')
    hold on
    %errorbar(readout_pos, conc_std_random(i, :), conc_SE_random(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'r')
    %errorbar(readout_pos, conc_std_average(i, :), conc_SE_average(i, :),'bo', 'LineWidth', LineWidth)      
    plot(readout_pos, conc_std_random(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'r')
    plot(readout_pos, conc_std_average(i, :), 'bo', 'LineWidth', LineWidth')
    xlabel(['Position x [µm]'])
    ylabel('\sigma_{C(x)}')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'linear')
    grid on
    
    % plot the relationship between CV_k and C_0
    subplot(3, 1, 3);
    errorbar(readout_pos, conc_CV_cilium(i, :), conc_CV_SE_cilium(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'g')
    hold on
    errorbar(readout_pos, conc_CV_random(i, :), conc_CV_SE_random(i, :), 'bo', 'LineWidth', LineWidth, 'Color', 'r')   
    errorbar(readout_pos, conc_CV_average(i, :), conc_CV_SE_average(i, :), 'bo', 'LineWidth', LineWidth)
    xlabel(['Position x [µm]'])
    ylabel(['CV_{C(x)} =  \sigma_{C(x)}/µ_{C(x)}'])
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'linear', 'YScale', 'log')
    grid on
   
end

 %% functions for the ODE
% reaction-diffusion equation
function dydx = odefun(x, y, c)
dC = -y(2,:) / D(c); % mass flux: j = -D*grad(C)
dj = p(c) * (c <= ncS) - d(c) * y(1,:); % conservation of mass: div(j) = p*H(-x) - d*C
dydx = [dC; dj];
end

% initial guess
function y = y0(x, c)
y = [0; 0];
end

% boundary & cell interface conditions
function res = bcfun(ya, yb)
res = ya(:);
res(1) = ya(2, 1); % zero flux at the left end of the source domain
res(2) = yb(2,nc); % zero flux at right end of the patterning domain
for c = 1:nc-1
    res(2*c+1) = ya(1,c+1) - yb(1,c); % concentration continuity
    res(2*c+2) = ya(2,c+1) - yb(2,c); % flux continuity
end
end

% log-normal distribution with adjusted mean & stddev
function pd = logndist(mu, sigma)
    pd = makedist('Lognormal', 'mu', log(mu/sqrt(1+(sigma/mu)^2)), 'sigma', sqrt(log(1+(sigma/mu)^2)));
end

end
