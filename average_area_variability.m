function average_area_variability 

% options
simulate = true; % if false, plot results from saved files instead of generating new data
write = true; % when creating new data, should the results be written to output files?
fitcosh = true; % fit an exp or cosh to the gradients
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 2; % number of independent simulation runs
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
%f8 = figure('Name', 'Dependency of gradient variability on cell size variability', 'Position', [0 0 2000 800]);

diameter = [0.49:0.5:4.9, 4.9:5:49.9]';

CV = 0.3; % fixed CV for {p, d, D, all}
CV_area = [0.5]; % mean cell area

%mu_a = (diameter/2)^2*pi;
%average_area_p = ncP * mu_a; % expected area if determinstic diameter for the patterning domain
%average_area_p
%average_area_s = ncS * mu_a; % expected area if determinstic diameter for the source domain
%plot_CV = CV_area(1);


% k = p, d, D, and all three together
% names = {'p', 'd', 'D', 'all'};
names = {'p'};
for k = 1:numel(names)
    
    filename = ['Cell_Area_vs_CV_' names{k} '.csv'];
    if names{k} == 'D'
        filename = 'Cell_Area_CV_Diff.csv';
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
        for i = 1:length(diameter)
            
            % calculate the area for varying diameters 
            mu_a = (diameter(i)/2)^2*pi;
            average_area_p = ncP * mu_a; % expected area patterning domain 
            average_area_s = ncS * mu_a; % expected area source domain 
            
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

                while sum(a_s_temp) < average_area_s 
                    
                    rand_area = random(logndist(mu_a, mu_a * CV_area), 1, 1);
                    
                    a_s_temp = [a_s_temp, rand_area];
                    
                end
                          
                assert(sum(a_s_temp(1:end-1)) < average_area_s, 'a_s_lower bigger than average area')
                assert(sum(a_s_temp) > average_area_s, 'a_s_temp smaller than average area')
                
                sum_s_upper = sum(a_s_temp);
                sum_s_lower = sum(a_s_temp(1:end-1));
                
                % edge case, where the source consists only of one cell 
                if length(a_s_temp) == 1
                
                    a_s = a_s_temp;
                    
                else    
                    if abs(sum_s_upper - average_area_s) < abs(sum_s_lower - average_area_s)
                        a_s = a_s_temp;
                    else
                        a_s = a_s_temp(1:end-1);
                    end 
                end
                
                              
                % from the area calculate the diameter of each cell 
                d_s = 2*sqrt((a_s)/pi);
                d_s_normalised = d_s / sum(d_s);
                d_s = d_s_normalised * LS;
                d_s = cumsum(d_s);                          
                d_s = fliplr(d_s);
                
                d_s

                % =========================================== %
                % Area computations for the Patterning Domain %
                % =========================================== %
                
                 % add cells as long as the overall arrea is not surpassed
                 while sum(a_p_temp) < average_area_p 
                    
                    rand_area = random(logndist(mu_a, mu_a * CV_area), 1, 1);
                    
                    a_p_temp = [a_p_temp, rand_area];
                    
                 end
                
                % add one more cell and check if the area is closer to the
                % desired area compaed to the case wihout this last cell
                 assert(sum(a_p_temp(1:end-1)) < average_area_p, 'a_p_lower bigger than average area')
                 assert(sum(a_p_temp) > average_area_p, 'a_p_temp smaller than average area')
                       
                sum_p_upper = sum(a_p_temp);
                sum_p_lower = sum(a_p_temp(1:end-1));
                        
                if abs(sum_p_upper - average_area_p) < abs(sum_p_lower - average_area_p)
                    a_p = a_p_temp;
                else
                    a_p = a_p_temp(1:end-1);
                end              
        
                % calculate the normalied diameter for each cell               
                d_p = 2*sqrt((a_p)/pi);
                d_p_normalised = d_p /sum(d_p);
                d_p = d_p_normalised * LP;
                d_p = cumsum(d_p);
                
            end

        end 
        
    end 
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
   
    
    




