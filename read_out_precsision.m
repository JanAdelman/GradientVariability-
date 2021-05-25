function readout_precision

% options
simulate = true; % if false, plot results from saved files instead of generating new data
write = true; % when creating new data, should the results be written to output files?
fitcosh = true; % fit an exp or cosh to the gradients
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 2; % number of independent simulation runs
nboot = 1e3; % number of bootstrap samples for error estimation
diameter = 4.9; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 19.26; % mean gradient length [µm]
mu_a = (diameter/2)^2*pi; %mean cell area 
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain
%CV = [0.01:0.01:0.09 0.1:0.05:0.95 1:0.5:10]'; % coefficient of variation of the kinetic parameters
CV = [0.0];
%plot_CV = CV(14); % plot the gradients for this CV value
CV_area = 0.5;
% analytical deterministic solution
nc = ncS + ncP; % total number of cells
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

CVfun = @(x) nanstd(x) / nanmean(x);
SEfun = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

fitopt = statset('TolFun', tol, 'TolX', tol);

close all
f1 = figure('Name', 'Individual gradients', 'Position', [0 0 2000 800]);
%names = {'p', 'd', 'D', 'all'};
names = {'p'};

% get the readout positions for the gradient 
readout_pos = [0:LP/10:LP];
%% vary molecular noise

f2 = figure('Name', 'Dependency of gradient parameters on molecular noise', 'Position', [0 0 2000 800]);
f3 = figure('Name', 'Dependency of gradient variability on molecular noise', 'Position', [0 0 2000 800]);
                    
% k = p, d, D, and all three together
for k = 1:numel(names)
    filename = ['CV_vs_CV_cv_area_0.5_' names{k} '.csv'];
    if names{k} == 'D'
        filename = 'CV_vs_CV_Diff_cv_area_0.5.csv';
    end
        
    if simulate
        lambda = NaN(length(CV), 1);
        lambda_SE = NaN(length(CV), 1);
        C0 = NaN(length(CV), 1);
        C0_SE = NaN(length(CV), 1);
        CV_lambda = NaN(length(CV), 1);
        CV_lambda_SE = NaN(length(CV), 1);
        CV_0 = NaN(length(CV), 1);
        CV_0_SE = NaN(length(CV), 1);

        % loop over variabilities
        for i = 1:length(CV)
            % loop over several independent runs
            fitted_lambda = NaN(nruns, 1);
            fitted_C0 = NaN(nruns, 1);
            
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
                 %assert(sum(l_p_temp(1:end-1)) < average_arel_p, 'l_p_lower bigger than average area')
                 %assert(sum(l_p_temp) > average_arel_p, 'l_p_temp smaller than average area')
                       
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
                
                % get initial solution 
                sol0 = bvpinit(x0, @y0);
                
                % solve the equation
                sol = bvp4c(@odefun, @bcfun, sol0, options);
                
                % initialise the start location at the beginning of ther
                % patterning domain 
                cell_beginning = 0;
                
                % array to store numerical integration results
                y_sol = [];
                % use the beginning of each cell as the x-coordinate 
                x_sol = [0, l_p(1:end-1)];
                              
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
                    
                    % get the average concentration per cell (thus
                    % normalised by cell length) 
                    trapz_sol = trapz(X,Y)/cell_length;
                    
                    % append the solution for each cell to the solution
                    % array 
                    y_sol = [y_sol, trapz_sol];
                    
                end 
                
                % fit an exponential in log space in the patterning domain
                param = polyfit(x_sol, log(y_sol), 1);
                fitted_lambda(j) = -1/param(1);
                fitted_C0(j) = exp(param(2));
                
                 fitted_lambda(j)
                 fitted_C0(j)

                % fit a hyperbolic cosine in log space in the patterning domain
                if fitcosh
                    logcosh = @(p,x) p(2) + log(cosh((LP-x)/p(1)));
                    mdl = fitnlm(x_sol, log(y_sol), logcosh, [fitted_lambda(j) log(fitted_C0(j)) - log(cosh(LP/fitted_lambda(j)))], 'Options', fitopt);
                    fitted_lambda(j) = mdl.Coefficients.Estimate(1);
                    fitted_C0(j) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP/fitted_lambda(j));
                    
                end
                
                % with the fitted values calculate the gradient and
                % evaluate, at different locations 
                c
                
                
                % plot the solution
                %{
                if CV(i) == plot_CV
                    figure(f1)
                    for s = 1:2
                        subplot(2, numel(names), k + (s-1)*numel(names))
                        hold all
                        plot(sol.x, sol.y(1,:), 'LineWidth', LineWidth)
                    end
                end
                %}
            end

            % determine the CV of the decay length and the amplitude over the independent runs
            % and also their standard errors from bootstrapping
            lambda(i) = mean(fitted_lambda);
            lambda_SE(i) = SEfun(fitted_lambda);
            C0(i) = mean(fitted_C0);
            C0_SE(i) = SEfun(fitted_C0);
            CV_lambda(i) = CVfun(fitted_lambda);
            CV_lambda_SE(i) = std(bootstrp(nboot, CVfun, fitted_lambda));
            CV_0(i) = CVfun(fitted_C0);
            CV_0_SE(i) = std(bootstrp(nboot, CVfun, fitted_C0));
        end
        
        % write data
        if write
            writetable(table(CV, lambda, lambda_SE, C0, C0_SE, CV_lambda, CV_lambda_SE, CV_0, CV_0_SE), filename);
        end
    else
        % read data
        T = readtable(filename);
        CV = T.CV;
        lambda = T.lambda;
        lambda_SE = T.lambda_SE;
        C0 = T.C0;
        C0_SE = T.C0_SE;
        CV_lambda = T.CV_lambda;
        CV_lambda_SE = T.CV_lambda_SE;
        CV_0 = T.CV_0;
        CV_0_SE = T.CV_0_SE;
    end

    % plot the relationship between CV_k and lambda
    figure(f2)
    subplot(2, numel(names), k)
    errorbar(CV, lambda, lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    xlabel(['CV_{' names{k} '}'])
    ylabel('\lambda [µm]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    grid on
    
    % plot the relationship between CV_k and C_0
    subplot(2, numel(names), k + numel(names))
    errorbar(CV, abs(C0-C(0)), C0_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    xlabel(['CV_{' names{k} '}'])
    ylabel('C_0 - \mu_0 [a.u.]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on
    
    % plot the relationship between CV_k and CV_lambda
    figure(f3)
    subplot(2, numel(names), k)
    errorbar(CV, CV_lambda, CV_lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    xlabel(['CV_{' names{k} '}'])
    ylabel('CV_\lambda')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on

    % plot the relationship between CV_k and CV_0
    subplot(2, numel(names), k + numel(names))
    errorbar(CV, CV_0, CV_0_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    xlabel(['CV_{' names{k} '}'])
    ylabel('CV_0')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on

    drawnow
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
