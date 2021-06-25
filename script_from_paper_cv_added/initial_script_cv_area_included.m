function gradient_variability_cv_area_adjusted

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
mu_a = (diameter/2)^2*pi; %mean cell area 
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain
CV = [0.01:0.01:0.09 0.1:0.05:0.95 1:0.5:10]'; % coefficient of variation of the kinetic parameters
plot_CV = CV(1); % plot the gradients for this CV value
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
names = {'p', 'd', 'D', 'all'};
%% vary molecular noise
%{
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
                
                % plot the solution
                if CV(i) == plot_CV
                    figure(f1)
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

%% vary domain length
ncP_arr = round(logspace(log10(20),log10(200),50))';
CV = 0.3;
f4 = figure('Name', 'Dependency of gradient parameters on domain size', 'Position', [0 0 2000 800]);
f5 = figure('Name', 'Dependency of gradient variability on domain size', 'Position', [0 0 2000 800]);
names = {'p', 'd', 'D', 'all'};
% k = p, d, D, and all three together
for k = 1:numel(names)
    
    filename = ['L_vs_CV_vs_cv_area_0.5_' names{k} '.csv'];
    if names{k} == 'D'
        filename = 'L_vs_cv_area_0.5_CV_Diff.csv';
    end
        
    if simulate
        L = ncP_arr * diameter;      
        lambda = NaN(length(ncP_arr), 1);
        lambda_SE = NaN(length(ncP_arr), 1);
        C0 = NaN(length(ncP_arr), 1);
        C0_SE = NaN(length(ncP_arr), 1);
        CV_lambda = NaN(length(ncP_arr), 1);
        CV_lambda_SE = NaN(length(ncP_arr), 1);
        CV_0 = NaN(length(ncP_arr), 1);
        CV_0_SE = NaN(length(ncP_arr), 1);

        % loop over patterning domain sizes
        for i = 1:length(ncP_arr)
            % solver initialization
            LP = L(i);
            
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
                
                % ======================================================= %
                % Solve the diffusion equation usiing  comstm grid vector
                % + update the cell counts for both comparments 
                % ======================================================= %
                                                           
                % create grid for the solver
                x0 = [];
                x0 = [x0, -l_s, 0, l_p];
   
                % initialise the solver
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
                
                nc = length(l_p) + length(l_s);
                ncS = length(l_s);
                ncP = length(l_p);
                               
                % default: all parameters constant
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
                
                options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);
                
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
                    LP = L(i);
                    logcosh = @(p,x) p(2) + log(cosh((LP-x)/p(1)));
                    mdl = fitnlm(sol.x(idx), log(sol.y(1,idx)), logcosh, [fitted_lambda(j) log(fitted_C0(j)) - log(cosh(LP/fitted_lambda(j)))], 'Options', fitopt);
                    fitted_lambda(j) = mdl.Coefficients.Estimate(1);
                    fitted_C0(j) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP/fitted_lambda(j));
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
        
        % write data
        if write
            writetable(table(L, lambda, lambda_SE, C0, C0_SE, CV_lambda, CV_lambda_SE, CV_0, CV_0_SE), filename);
        end
    else
        % read data
        T = readtable(filename);
        L = T.L;
        lambda = T.lambda;
        lambda_SE = T.lambda_SE;
        C0 = T.C0;
        C0_SE = T.C0_SE;
        CV_lambda = T.CV_lambda;
        CV_lambda_SE = T.CV_lambda_SE;
        CV_0 = T.CV_0;
        CV_0_SE = T.CV_0_SE;
    end
    
    % plot the relationship between L and lambda    
    figure(f4)
    subplot(2, numel(names), k)
    errorbar(L, lambda, lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title(['CV_{' names{k} '} = ' num2str(CV) ', CV_area = ' num2str(CV_area)])
    xlabel('L [µm]')
    ylabel('\lambda [µm]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    
    % plot the relationship between L and C_0
    subplot(2, numel(names), k + numel(names))
    errorbar(L, C0, C0_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title(['CV_{' names{k} '} = ' num2str(CV) ', CV_area = ' num2str(CV_area)])
    xlabel('L [µm]')
    ylabel('C_0 [a.u.]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    
    % plot the relationship between L and CV_lambda
    figure(f5)
    subplot(2, numel(names), k)
    errorbar(L, CV_lambda, CV_lambda_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title(['CV_{' names{k} '} = ' num2str(CV) ', CV_area = ' num2str(CV_area)])
    xlabel('L [µm]')
    ylabel('CV_\lambda')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    grid on

    % plot the relationship between L and CV_0
    subplot(2, numel(names), k + numel(names))
    errorbar(L, CV_0, CV_0_SE, 'bo', 'LineWidth', LineWidth)
    hold on
    title(['CV_{' names{k} '} = ' num2str(CV) ', CV_area = ' num2str(CV_area)])
    xlabel('L [µm]')
    ylabel('CV_0')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    
    drawnow
end
%}
%% vary white noise
CV_white = logspace(-8,0,161)'; % relative strength of white noise

plot_CV = 0.1;
CV = 0.3;
CV_area = 0.5;

% solver initialization
ncP = 50;
nc = ncS + ncP;
LP = ncP * diameter;
           
f6 = figure('Name', 'Dependency of gradient parameters on white noise', 'Position', [0 0 2000 800]);
f7 = figure('Name', 'Dependency of gradient variability on white noise', 'Position', [0 0 2000 800]);
x = (-ncS+0.5:ncP-0.5) * diameter;


%x = (-LS+0.5:4.9:(LP-0.5))
% k = p, d, D, and all three together

%names = {'p', 'd', 'D', 'all'};
names = {'all'};

for k = 1:numel(names)
    
    filename = ['SNR_vs_CV_vs_CV_area_' names{k} '.csv'];
    
    if names{k} == 'D'
        filename = 'script_from_paper_cv_added/SNR_vs_CV_vs_CV_area_Diff.csv';

    end
    
    if simulate
        lambda = NaN(length(CV_white), 1);
        lambda_SE = NaN(length(CV_white), 1);
        C0 = NaN(length(CV_white), 1);
        C0_SE = NaN(length(CV_white), 1);
        CV_lambda = NaN(length(CV_white), 1);
        CV_lambda_SE = NaN(length(CV_white), 1);
        CV_0 = NaN(length(CV_white), 1);
        CV_0_SE = NaN(length(CV_white), 1);
        
        lambda_lin = lambda;
        lambda_SE_lin = lambda_SE;
        C0_lin = C0;
        C0_SE_lin = C0_SE;
        CV_lambda_lin = CV_lambda;
        CV_lambda_SE_lin = CV_lambda_SE;
        CV_0_lin = CV_0;
        CV_0_SE_lin = CV_0_SE;

        % loop over patterning domain sizes
        for i = 1:length(CV_white)
            % loop over several independent runs
            fitted_lambda = NaN(nruns, 1);
            fitted_C0 = NaN(nruns, 1);
            fitted_lambda_lin = fitted_lambda;
            fitted_C0_lin = fitted_C0;
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
                              
                % default: all parameters constant
               
                p = mu_p * ones(nc, 1);
                d = mu_d * ones(nc, 1);
                D = mu_D * ones(nc, 1);

                % draw random kinetic parameters for each cell
                %{
                if k == 1 || k == 4
                    p = random(logndist(mu_p, mu_p * CV), nc, 1);
                end
                if k == 2 || k == 4
                    d = random(logndist(mu_d, mu_d * CV), nc, 1);
                end
                if k == 3 || k == 4
                    D = random(logndist(mu_D, mu_D * CV), nc, 1);
                end
                %}
                p = random(logndist(mu_p, mu_p * CV), nc, 1);
                d = random(logndist(mu_d, mu_d * CV), nc, 1);
                D = random(logndist(mu_D, mu_D * CV), nc, 1);
                
                % solve the equation
                options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);
                
                 % get initial solution 
                sol0 = bvpinit(x0, @y0);
                
                sol = bvp4c(@odefun, @bcfun, sol0, options);
                             
                % add white noise
                y = deval(sol, x);
             
                %y = y(1,:) + normrnd(0, C(0) * CV_white(i), [1 nc]);
                y = y(1,:) + normrnd(0, C(0) * CV_white(i), [1 length(y)]);
              
                % fit an exponential in log space in the patterning domain
                idx = find(x >= 0 & y > 0);
                param = polyfit(x(idx), log(y(idx)), 1);
                fitted_lambda(j) = -1/param(1);
                fitted_C0(j) = exp(param(2));
                
                % fit an exponential in lin space in the patterning domain
                try
                    idx = find(x >= 0);
                    linexp = @(p,x) exp(polyval(p,x));
                    mdl = fitnlm(x(idx), y(idx), linexp, [-1/fitted_lambda(j) 0], 'Options', fitopt);
                    fitted_lambda_lin(j) = -1/mdl.Coefficients.Estimate(1);
                    fitted_C0_lin(j) = exp(mdl.Coefficients.Estimate(2));
                end
                
                % fit a hyperbolic cosine in log space in the patterning domain
                if fitcosh
                    try
                        idx = find(x >= 0 & y > 0);
                        logcosh = @(p,x) p(2) + log(cosh((LP-x)/p(1)));
                        mdl = fitnlm(x(idx), log(y(idx)), logcosh, [fitted_lambda(j) log(fitted_C0(j)) - log(cosh(LP/fitted_lambda(j)))], 'Options', fitopt);
                        fitted_lambda(j) = mdl.Coefficients.Estimate(1);
                        fitted_C0(j) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP/fitted_lambda(j));
                    end
                end
                
                % plot the solution
                if CV_white(i) == plot_CV
                    figure(f1)
                    for s = 1:2
                        subplot(2, numel(names), k + (s-1)*numel(names))
                        hold all
                        plot(x, y(1,:), 'LineWidth', LineWidth)
                    end
                end
                
            end

            % determine the CV of the decay length and the amplitude over the independent runs
            % and also their standard errors from bootstrapping
            lambda(i) = nanmean(fitted_lambda);
            lambda_SE(i) = SEfun(fitted_lambda);
            C0(i) = nanmean(fitted_C0);
            C0_SE(i) = SEfun(fitted_C0);
            CV_lambda(i) = CVfun(fitted_lambda);
            CV_lambda_SE(i) = nanstd(bootstrp(nboot, CVfun, fitted_lambda));
            CV_0(i) = CVfun(fitted_C0);
            CV_0_SE(i) = nanstd(bootstrp(nboot, CVfun, fitted_C0));
            
            lambda_lin(i) = nanmean(fitted_lambda_lin);
            lambda_SE_lin(i) = SEfun(fitted_lambda_lin);
            C0_lin(i) = nanmean(fitted_C0_lin);
            C0_SE_lin(i) = SEfun(fitted_C0_lin);
            CV_lambda_lin(i) = CVfun(fitted_lambda_lin);
            CV_lambda_SE_lin(i) = nanstd(bootstrp(nboot, CVfun, fitted_lambda_lin));
            CV_0_lin(i) = CVfun(fitted_C0_lin);
            CV_0_SE_lin(i) = nanstd(bootstrp(nboot, CVfun, fitted_C0_lin));
        end
 
        % write data
        if write
            writetable(table(CV_white, lambda, lambda_SE, C0, C0_SE, CV_lambda, CV_lambda_SE, CV_0, CV_0_SE, lambda_lin, lambda_SE_lin, C0_lin, C0_SE_lin, CV_lambda_lin, CV_lambda_SE_lin, CV_0_lin, CV_0_SE_lin), filename);
        end
    else
        % read data
        T = readtable(filename);
        CV_white = T.CV_white;
        lambda = T.lambda;
        lambda_SE = T.lambda_SE;
        C0 = T.C0;
        C0_SE = T.C0_SE;
        CV_lambda = T.CV_lambda;
        CV_lambda_SE = T.CV_lambda_SE;
        CV_0 = T.CV_0;
        CV_0_SE = T.CV_0_SE;
        
        lambda_lin = T.lambda_lin;
        lambda_SE_lin = T.lambda_SE_lin;
        C0_lin = T.C0_lin;
        C0_SE_lin = T.C0_SE_lin;
        CV_lambda_lin = T.CV_lambda_lin;
        CV_lambda_SE_lin = T.CV_lambda_SE_lin;
        CV_0_lin = T.CV_0_lin;
        CV_0_SE_lin = T.CV_0_SE_lin;
    end
    
    % plot the relationship between CV_white and lambda
    figure(f6)
    subplot(2, numel(names), k)
    hold on
    box on
    grid on
    errorbar(CV_white, lambda, lambda_SE, 'bo', 'LineWidth', LineWidth)
    errorbar(CV_white, lambda_lin, lambda_SE_lin, 'ro', 'LineWidth', LineWidth)
    title(['CV_{' names{k} '} = ' num2str(CV)])
    xlabel('CV_{white}')
    ylabel('\lambda [µm]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    xticks(10.^(-8:2:0))
    xlim([min(CV_white) max(CV_white)])
    ylim([0 100])
    
    % plot the relationship between CV_white and C_0
    subplot(2, numel(names), k + numel(names))
    hold on
    box on
    grid on
    errorbar(CV_white, C0, C0_SE, 'bo', 'LineWidth', LineWidth)
    errorbar(CV_white, C0_lin, C0_SE_lin, 'ro', 'LineWidth', LineWidth)
    title(['CV_{' names{k} '} = ' num2str(CV)])
    xlabel('CV_{white}')
    ylabel('C_0 [a.u.]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    xticks(10.^(-8:2:0))
    xlim([min(CV_white) max(CV_white)])
    ylim([0 0.6])
    
    drawnow
    
    % plot the relationship between CV_white and CV_lambda
    figure(f7)
    subplot(2, numel(names), k)
    hold on
    box on
    grid on
    errorbar(CV_white, CV_lambda, CV_lambda_SE, 'bo', 'LineWidth', LineWidth)
    errorbar(CV_white, CV_lambda_lin, CV_lambda_SE_lin, 'ro', 'LineWidth', LineWidth)
    title(['CV_{' names{k} '} = ' num2str(CV)])
    xlabel('CV_{white}')
    ylabel('CV_\lambda')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
    xticks(10.^(-8:2:0))
    xlim([min(CV_white) max(CV_white)])
    ylim([0.01 1])
    
    % plot the relationship between CV_white and CV_0
    subplot(2, numel(names), k + numel(names))
    hold on
    box on
    grid on
    errorbar(CV_white, CV_0, CV_0_SE, 'bo', 'LineWidth', LineWidth)
    errorbar(CV_white, CV_0_lin, CV_0_SE_lin, 'ro', 'LineWidth', LineWidth)
    title(['CV_{' names{k} '} = ' num2str(CV)])
    xlabel('CV_{white}')
    ylabel('CV_0')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log')
    xticks(10.^(-8:2:0))
    xlim([min(CV_white) max(CV_white)])
    ylim([0 0.5])
    
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

