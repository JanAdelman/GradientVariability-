function readout_precision
% options
simulate = true; % if false, plot results from saved files instead of generating new data
write = true; % when creating new data, should the results be written to output files?
fitcosh = true; % fit an exp or cosh to the gradients
LineWidth = 100;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 4.9; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 19.26; % mean gradient length [µm]
mu_a = (diameter/2)^2*pi; %mean cell area 
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain
CV = [0.3]; % coefficient of variation of the kinetic parameters
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

hill = @(x,k_hill) 1*(x ./ (x + k_hill));

% high resolution k values
k_d = logspace(log10(C(0)), log10(C(LP)), 100);

CV_K = [0, 0.001, 0.1, 0.3];

filename_pos_neg_average_lin = ['CV_K_Parameter_Tables_Hill/pos_neg_average_lin_hill.csv'];
filename_pos_neg_average = ['CV_K_Parameter_Tables_Hill/pos_neg_average_hill.csv'];
filename_pos_neg_cilium_lin = ['CV_K_Parameter_Tables_Hill/pos_neg_cilium_lin_hill.csv'];
filename_pos_neg_cilium = ['CV_K_Parameter_Tables_Hill/pos_neg_cilium_hill.csv'];
filename_pos_neg_random_lin = ['CV_K_Parameter_Tables_Hill/pos_neg_random_lin_hill.csv'];
filename_pos_neg_random = ['CV_K_Parameter_Tables_Hill/pos_neg_random_hill.csv'];

pos_neg_average_lin = NaN(length(CV_K), 2);
pos_neg_average = NaN(length(CV_K), 2);
pos_neg_cilium_lin = NaN(length(CV_K), 2);
pos_neg_cilium = NaN(length(CV_K), 2);
pos_neg_random_lin = NaN(length(CV_K), 2);
pos_neg_random = NaN(length(CV_K), 2);



for k = 1:numel(CV_K)
    
        ave_arr_lin_low = ones(length(k_d), 1);
        ave_arr_low = ones(length(k_d), 1);
        cil_arr_lin_low = ones(length(k_d), 1);
        cil_arr_low = ones(length(k_d), 1);
        rand_arr_low = ones(length(k_d), 1);
        rand_arr_lin_low = ones(length(k_d), 1);
    
        ave_arr_lin_high = ones(length(k_d), 1);
        ave_arr_high = ones(length(k_d), 1);
        cil_arr_lin_high = ones(length(k_d), 1);
        cil_arr_high = ones(length(k_d), 1);
        rand_arr_high = ones(length(k_d), 1);
        rand_arr_lin_high = ones(length(k_d), 1);
   
        filename_readout_pos_average = ['CV_K_Parameter_Tables_Hill/readout_pos_average_' num2str(CV_K(k)) '_hill.csv'];
        filename_readout_pos_cilium = ['CV_K_Parameter_Tables_Hill/readout_pos_cilium_' num2str(CV_K(k)) '_hill.csv'];
        filename_readout_pos_random = ['CV_K_Parameter_Tables_Hill/readout_pos_random_' num2str(CV_K(k)) '_hill.csv'];
        filename_readout_pos_average_lin = ['CV_K_Parameter_Tables_Hill/readout_pos_average_' num2str(CV_K(k)) '_hill_lin_inter.csv'];
        filename_readout_pos_cilium_lin = ['CV_K_Parameter_Tables_Hill/readout_pos_cilium_' num2str(CV_K(k)) '_hill_lin_inter.csv'];
        filename_readout_pos_random_lin = ['CV_K_Parameter_Tables_Hill/readout_pos_random_' num2str(CV_K(k)) '_hill_lin_inter.csv'];     

        % allocatte memory to store the interp. values 
        half_max_average_array = NaN(nruns, length(k_d));
        half_max_cilium_array = NaN(nruns, length(k_d));
        half_max_random_array = NaN(nruns, length(k_d));
        half_max_average_array_lin = NaN(nruns, length(k_d));
        half_max_cilium_array_lin = NaN(nruns, length(k_d));
        half_max_random_array_lin = NaN(nruns, length(k_d));
               
        % loop over variabilities
        for i = 1:length(k_d)
            
            half_max_average_array_one_nun = NaN(1, nruns);
            half_max_cilium_array_one_run = NaN(1, nruns);
            half_max_random_array_one_run = NaN(1, nruns);
            half_max_average_array_lin_one_run = NaN(1, nruns);
            half_max_cilium_array_lin_one_run = NaN(1, nruns);
            half_max_random_array_lin_one_run  = NaN(1, nruns);
            
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

                p = random(logndist(mu_p, mu_p * CV(1)), nc, 1);
                d = random(logndist(mu_d, mu_d * CV(1)), nc, 1); 
                D = random(logndist(mu_D, mu_D * CV(1)), nc, 1);
                
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
                
                rand_x_points =[];
                                
                for cell_loc = 1:ncP
                    
                    % set the upper interval as the end of a cell 
                    cell_end = l_p(cell_loc);
                    
                    % define interval where to extract solutions 
                    logical_indexes = (sol.x <= cell_end) & (sol.x >= cell_beginning);
                    % extract indices of the desired solutions 
                    interval = find(logical_indexes);
                    
                    % get lenght of the cell for normalisation
                    cell_length = cell_end - cell_beginning;
                    
                    % get the x and y solution 
                    X = sol.x(interval);
                    Y = sol.y(1, interval);
                    
                    % get unique x, y values for interpolation solver
                    x_unique = unique(X,'stable');
                    y_unique = unique(Y,'stable');
                    
                    % closest to mean 
                    %[min_distance, index] = min( abs(X - mid_point(cell_loc)));
                    %cilium_sol = Y(index);
                    
                    % interpolation um Mittelwert zu bestimmen 
                    cilium_sol = interp1(x_unique, y_unique, mid_point(cell_loc));
                    
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_cilium = [y_sol_cilium, cilium_sol];
                    
                    % get a concentration randomly from a cell-point 
                    % and use this point with interpolation 
                    %rand_conc = randsample(Y,1);
                    
                    % get a random point in the cell
                    rand_x = (cell_end-cell_beginning).*rand(1,1) + cell_beginning;
                    
                    rand_x_points = [rand_x_points, rand_x];
                                       
                    % get solution at that point
                    rand_conc = interp1(x_unique, y_unique, rand_x);
                    
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_random = [y_sol_random, rand_conc];
                    
                    % get the average concentration per cell (thus
                    % normalised by cell length) 
                    trapz_sol = trapz(X,Y)/cell_length;
  
                    % append the solution for each cell to the solution
                    % array 
                    y_sol_average = [y_sol_average, trapz_sol];
                    
                    % set the lower interval for the next iteration as the
                    % current end of the cell 
                    cell_beginning = cell_end;
                    
                                      
                end 
              
                k_d_iter = k_d(i);
                
                % add noise to system 
                k_noisy = random(logndist(k_d_iter, k_d_iter * CV_K(k)), 1, 1);
                
                hill_random_point = hill(y_sol_random, k_noisy);
                hill_average_conc = hill(y_sol_average, k_noisy);
                hill_cilium = hill(y_sol_cilium, k_noisy);
           
                half_max_average_lin = interp1(hill_average_conc, x_sol, 0.5, 'linear','extrap');
                half_max_random_lin = interp1(hill_random_point, x_sol, 0.5, 'linear','extrap');
                half_max_cilium_lin = interp1(hill_cilium, x_sol, 0.5, 'linear','extrap');
                           
                % cubic interplation 
                half_max_average = pchip(hill_average_conc, x_sol, 0.5);
                half_max_random = pchip(hill_random_point, rand_x_points, 0.5);
                half_max_cilium = pchip(hill_cilium, mid_point, 0.5);
                
                if half_max_average_lin < 0
                    ave_arr_lin_low(i, 1) = 2;
                end
                if half_max_cilium_lin < 0
                    cil_arr_lin_low(i, 1) = 2;
                end
                if half_max_random_lin < 0
                    rand_arr_lin_low(i, 1) = 2;
                end
                if half_max_average < 0
                    ave_arr_low(i, 1) = 2;
                end
                if half_max_cilium < 0
                    cil_arr_low(i, 1) = 2;
                end
                if half_max_random < 0
                    rand_arr_low(i, 1) = 2;
                end
                
                if ((half_max_average > 245))
                    ave_arr_high(i, 1) = 2;
                    
                end
                
                if ((half_max_cilium > 245) )
                    cil_arr_high(i, 1) = 2;
        
                end
                
                if ((half_max_random > 245))
                    rand_arr_high(i, 1) = 2;
               
                end
                
                if ((half_max_average_lin > 245))
                    ave_arr_lin_high(i, 1) = 2;
              
                end
                
                if ((half_max_cilium_lin > 245))
                    cil_arr_lin_high(i, 1) = 2;
                    
                end
                
                if ((half_max_random_lin > 245))
                    rand_arr_lin_high(i, 1) = 2;
                
                end
               
                               
                half_max_average_array_one_run(1,j)  = half_max_average;
                half_max_cilium_array_one_run(1,j) = half_max_cilium;
                half_max_random_array_one_run(1,j) = half_max_random;
                half_max_average_array_lin_one_run(1,j) = half_max_average_lin;
                half_max_cilium_array_lin_one_run(1,j) = half_max_cilium_lin;
                half_max_random_array_lin_one_run(1,j)  = half_max_random_lin;
                              
                
            end
           
            half_max_average_array(:, i) =  half_max_average_array_one_run';
            half_max_cilium_array(:, i) = half_max_cilium_array_one_run';
            half_max_random_array(:, i) = half_max_random_array_one_run';
            half_max_average_array_lin(:, i) =  half_max_average_array_lin_one_run';
            half_max_cilium_array_lin(:, i) = half_max_cilium_array_lin_one_run';
            half_max_random_array_lin(:, i) = half_max_random_array_lin_one_run';
              
           
        end
        
        mean_x_average = nanmean(half_max_average_array);
        std_x_average = nanstd(half_max_average_array);
        SE_x_average = nanstd(bootstrp(nboot, SEfun, half_max_average_array));
           
        mean_x_cilium = nanmean(half_max_cilium_array); 
        std_x_cilium = nanstd(half_max_cilium_array);
        SE_x_cilium = nanstd(bootstrp(nboot, SEfun, half_max_cilium_array));
           
        mean_x_random = nanmean(half_max_random_array);
        std_x_random = nanstd(half_max_random_array);
        SE_x_random = nanstd(bootstrp(nboot, SEfun, half_max_random_array));
           
        mean_x_average_lin = nanmean(half_max_average_array_lin);
        std_x_average_lin = nanstd(half_max_average_array_lin);
        SE_x_average_lin = nanstd(bootstrp(nboot, SEfun, half_max_average_array_lin));
                   
        mean_x_cilium_lin = nanmean(half_max_cilium_array_lin); 
        std_x_cilium_lin = nanstd(half_max_cilium_array_lin);
        SE_x_cilium_lin = nanstd(bootstrp(nboot, SEfun, half_max_cilium_array_lin));
           
        mean_x_random_lin = nanmean(half_max_random_array_lin);
        std_x_random_lin = nanstd(half_max_random_array_lin);
        SE_x_random_lin = nanstd(bootstrp(nboot, SEfun, half_max_random_array_lin));
        
               
        % write data
        if write
             writetable(table(k_d', mean_x_average', std_x_average',  SE_x_average', 'VariableNames', {'k_d', 'mean_pos_average', 'std_pos_average', 'SE_std'}), filename_readout_pos_average);
             writetable(table(k_d', mean_x_random', std_x_random',SE_x_random', 'VariableNames', {'k_d', 'mean_pos_random', 'std_pos_random', 'SE_std'}), filename_readout_pos_random);
             writetable(table(k_d', mean_x_cilium', std_x_cilium',SE_x_cilium', 'VariableNames', {'k_d', 'mean_pos_cilium', 'std_pos_cilium', 'SE_std'}), filename_readout_pos_cilium);
             writetable(table(k_d', mean_x_average_lin', std_x_average_lin', SE_x_average_lin', 'VariableNames', {'k_d', 'mean_pos_average', 'std_pos_average', 'SE_std'}), filename_readout_pos_average_lin);
             writetable(table(k_d', mean_x_random_lin', std_x_random_lin', SE_x_random_lin', 'VariableNames', {'k_d', 'mean_pos_random', 'std_pos_random', 'SE_std'}), filename_readout_pos_random_lin);
             writetable(table(k_d', mean_x_cilium_lin', std_x_cilium_lin',SE_x_cilium_lin',  'VariableNames', {'k_d', 'mean_pos_cilium', 'std_pos_cilium', 'SE_std'}), filename_readout_pos_cilium_lin);
        end
  
        index_average_lin = max(find(ave_arr_lin_low == 2));
        index_average =  max(find(ave_arr_low == 2));
        index_cilium_lin = max(find(cil_arr_lin_low == 2));
        index_cilium = max(find(cil_arr_low  == 2));
        index_random_lin = max(find(rand_arr_lin_low  == 2));
        index_random = max(find(rand_arr_low  == 2));
    
        min_x_average_lin = mean_x_average_lin(index_average_lin);
        min_x_average = mean_x_average(index_average);
        min_x_cilium_lin = mean_x_cilium_lin(index_cilium_lin);
        min_x_cilium = mean_x_cilium(index_cilium);
        min_x_random_lin = mean_x_random_lin(index_random_lin);
        min_x_random = mean_x_random(index_random);
    
     index_average_lin_high = min(find(ave_arr_lin_high  == 2 ));
     index_average_high =  min(find(ave_arr_high  == 2 ));
     index_cilium_lin_high = min(find(cil_arr_lin_high == 2));
     index_cilium_high = min(find(cil_arr_high == 2));
     index_random_lin_high = max(find(rand_arr_lin_high == 2));
     index_random_high = min(find(rand_arr_high == 2));
        
        
     max_x_average_lin = mean_x_average_lin(index_average_lin_high);
     max_x_average = mean_x_average(index_average_high);
     max_x_cilium_lin = mean_x_cilium_lin(index_cilium_lin_high);
     max_x_cilium = mean_x_cilium(index_cilium_high);
     max_x_random_lin = mean_x_random_lin(index_random_lin_high);
     max_x_random = mean_x_random(index_random_high);
    
     pos_neg_average_lin(k, :)  =  [min_x_average_lin, max_x_average_lin];
     pos_neg_average(k, :) = [min_x_average, max_x_average];
     pos_neg_cilium_lin(k, :)  =  [min_x_cilium_lin,  max_x_cilium_lin];
     pos_neg_cilium(k, :) = [min_x_cilium,  max_x_cilium];
     pos_neg_random_lin(k, :)  =  [min_x_random_lin, max_x_random_lin];
     pos_neg_random(k, :) = [min_x_random, max_x_random];
    end

writetable(table(CV_K', pos_neg_average_lin(:,1), pos_neg_average_lin(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_average_lin);
writetable(table(CV_K', pos_neg_average(:,1), pos_neg_average(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_average);
writetable(table(CV_K', pos_neg_cilium_lin(:,1), pos_neg_cilium_lin(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_cilium_lin);
writetable(table(CV_K', pos_neg_cilium(:,1), pos_neg_cilium(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_cilium);
writetable(table(CV_K', pos_neg_random_lin(:,1), pos_neg_random_lin(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_random_lin);
writetable(table(CV_K', pos_neg_random(:,1), pos_neg_random(:,2), 'VariableNames', {'CV_K', 'min_pos', 'max_pos'}), filename_pos_neg_random);

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
