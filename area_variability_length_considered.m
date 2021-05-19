function area_variability_2
% set seed
%rng('default');
%rng(10);

% options
simulate = true; % if false, plot results from saved files instead of generating new data
write = true; % when creating new data, should the results be written to output files?
fitcosh = true; % fit an exp or cosh to the gradients
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 4.9; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 19.26; % mean gradient length [µm]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain

% analytical deterministic solution
nc = ncS + ncP; % total number of cells
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

CVfun = @(x) nanstd(x) / nanmean(x);
SEfun = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

fitopt = statset('TolFun', tol, 'TolX', tol);

close all
f10 = figure('Name', 'Individual gradients', 'Position', [0 0 2000 800]);

%% vary cell size 
f7 = figure('Name', 'Dependency of gradient parameters on cell size variability', 'Position', [0 0 2000 800]);
f8 = figure('Name', 'Dependency of gradient variability on cell size variability', 'Position', [0 0 2000 800]);


CV = 0.3; % fixed CV for {p, d, D, all}
CV_area = [0.01:0.01:0.09 0.1:0.05:0.95 1:0.5:10]';
mu_a = (diameter/2)^2*pi; % mean cell area
average_area_p = ncP * mu_a; % expected area if determinstic diameter for the patterning domain
average_area_s = ncS * mu_a; % expected area if determinstic diameter for the source domain
%plot_CV = CV_area(14);

plot_CV = CV_area(1);
average_length_p = ncP * diameter;
average_length_s = ncS * diameter;

% k = p, d, D, and all three together
names = {'p', 'd', 'D', 'all'};
%names = {'p'}
for k = 1:numel(names)
    
    filename = ['Cell_Area_vs_CV_length_version_' names{k} '.csv'];
    if names{k} == 'D'
        filename = 'Cell_Area_CV_Diff_length_version.csv';
    end
    
    if simulate
        
        lambda = NaN(length(CV_area), 1);
        lambda_SE = NaN(length(CV_area), 1);
        C0 = NaN(length(CV_area), 1);
        C0_SE = NaN(length(CV_area), 1);
        CV_lambda = NaN(length(CV_area), 1);
        CV_lambda_SE = NaN(length(CV_area), 1);
        CV_0 = NaN(length(CV_area), 1);
        CV_0_SE = NaN(length(CV_area), 1);

        % loop over patterning domain sizes
        for i = 1:length(CV_area)
            
            fitted_lambda = NaN(nruns, 1);
            fitted_C0 = NaN(nruns, 1);
            
            % loop over several independent runs
            for j = 1:nruns
               
                % initialise arrays for variable cell size          
                a_s_temp = [];
                a_s = [];
                            
                a_p_temp = [];
                a_p = [];
                
                % ======================================= %
                % Area computations for the Source Domain %
                % ======================================= %

                while sum(a_s_temp) < average_length_s 
                    
                    rand_area = random(logndist(mu_a, mu_a * CV_area(i)), 1, 1);
                    length_s_sampled = 2*sqrt((rand_area)/pi);
                    a_s_temp = [a_s_temp, length_s_sampled];
                                        
                end
                       
                %assert(sum(a_s_temp(1:end-1)) < average_area_s, 'a_s_lower bigger than average area')
                %assert(sum(a_s_temp) > average_area_s, 'a_s_temp smaller than average area')
                
                sum_s_upper = sum(a_s_temp);
                sum_s_lower = sum(a_s_temp(1:end-1));
                
                % edge case, where the source consists only of one cell 
                if length(a_s_temp) == 1
                
                    a_s = a_s_temp;
                    
                else    
                    if abs(sum_s_upper - average_length_s) < abs(sum_s_lower - average_length_s)
                        a_s = a_s_temp;
                    else
                        a_s = a_s_temp(1:end-1);
                    end 
                end
                
                              
                % from the area calculate the diameter of each cell 
                d_s_normalised = a_s / sum(a_s);
                a_s = d_s_normalised * LS;
                a_s = cumsum(a_s);                          
                a_s = fliplr(a_s);
                

                % =========================================== %
                % Area computations for the Patterning Domain %
                % =========================================== %
                
                 % add cells as long as the overall arrea is not surpassed
                 while sum(a_p_temp) < average_length_p 

                    rand_area = random(logndist(mu_a, mu_a * CV_area(i)), 1, 1);
                    
                    length_p_sampled = 2*sqrt((rand_area)/pi);
                                      
                    a_p_temp = [a_p_temp, length_p_sampled];
                    
                 end
                
                % add one more cell and check if the area is closer to the
                % desired area compaed to the case wihout this last cell
                 %assert(sum(a_p_temp(1:end-1)) < average_area_p, 'a_p_lower bigger than average area')
                 %assert(sum(a_p_temp) > average_area_p, 'a_p_temp smaller than average area')
                       
                sum_p_upper = sum(a_p_temp);
                sum_p_lower = sum(a_p_temp(1:end-1));
                        
                if abs(sum_p_upper - average_length_p) < abs(sum_p_lower - average_length_p)
                    a_p = a_p_temp;
                else
                    a_p = a_p_temp(1:end-1);
                end              
        
                % calculate the normalied diameter for each cell               
                
                d_p_normalised = a_p /sum(a_p);
                a_p = d_p_normalised * LP;
                a_p = cumsum(a_p);
                
                
                
                % ======================================================= %
                % Solve the diffusion equation usiing  comstm grid vector
                % + update the cell counts for both comparments 
                % ======================================================= %
                
                                           
                % create grid for the solver
                x0 = [];
                x0 = [x0, -a_s, 0, a_p];
                
                % initialise the solver
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
                
                nc = length(a_p) + length(a_s);
                ncS = length(a_s);
                ncP = length(a_p);

               
                options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);
                
                % allocate memory fr the kineetic parameters 
                p = mu_p * ones(nc, 1);
                d = mu_d * ones(nc, 1);
                D = mu_D * ones(nc, 1);
                
                
                % draw random kinetic parameters for each cell
                if k == 1 || k == 4
                    p = random(logndist(mu_p, mu_p * CV), nc, 1);
                end
                if k == 2 || k == 4
                    d = random(logndist(mu_d, mu_d * CV), nc, 1);
                end
                if k == 3 || k == 4
                    D = random(logndist(mu_D, mu_D * CV), nc, 1);
                end
                
                % get initial solution
                sol0 = bvpinit(x0, @y0);
                
                % solve the equation
                sol = bvp4c(@odefun, @bcfun, sol0, options);

                % fit an exponential in log space in the patterning domain
                idx = find(sol.x >= 0);
                
                param = polyfit(sol.x(idx), log(sol.y(1,idx)), 1);
                fitted_lambda(j) = -1/param(1);
                fitted_C0(j) = exp(param(2));
                   
                % fit a hyperbolic cosine in log space in the patterning domain
                
                if fitcosh
                    logcosh = @(p,x) p(2) + log(cosh((LP-x)/p(1)));
                    mdl = fitnlm(sol.x(idx), log(sol.y(1,idx)), logcosh, [fitted_lambda(j) log(fitted_C0(j)) - log(cosh(LP/fitted_lambda(j)))], 'Options', fitopt);
                    fitted_lambda(j) = mdl.Coefficients.Estimate(1);
                    fitted_C0(j) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP/fitted_lambda(j));
                end
                
                if CV_area(i) == plot_CV
                    figure(f10)
                    for s = 1:2
                        subplot(2, numel(names), k + (s-1)*numel(names))
                        hold all
                        plot(sol.x, sol.y(1,:), 'LineWidth', LineWidth)
                    end
                end
                
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
        
        if write
            writetable(table(CV_area, lambda, lambda_SE, C0, C0_SE, CV_lambda, CV_lambda_SE, CV_0, CV_0_SE), filename);
        end        
        
    end
    
    % plot the relationship between CV_area and lambda
    figure(f7)
    subplot(2, numel(names), k)
    errorbar(CV_area, lambda, lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title((['CV_{' names{k} '}']))
    xlabel(['CV_{area}'])
    ylabel('\lambda [µm]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    grid on
    
    % plot the relationship between CV_area and C_0
    subplot(2, numel(names), k + numel(names))
    errorbar(CV_area, abs(C0-C(0)), C0_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    xlabel(['CV_{area}'])
    ylabel('C_0 - \mu_0 [a.u.]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on
    
     % plot the relationship between CV_k and CV_lambda
    figure(f8)
    subplot(2, numel(names), k)
    errorbar(CV_area, CV_lambda, CV_lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title((['CV_{' names{k} '}']))
    xlabel(['CV_{area}'])
    ylabel('CV_\lambda')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on

    % plot the relationship between CV_k and CV_0
    subplot(2, numel(names), k + numel(names))
    errorbar(CV_area, CV_0, CV_0_SE, 'bo', 'LineWidth', LineWidth)
    hold on   
    xlabel(['CV_{area}'])
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
