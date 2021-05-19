function test

ncS = 5; % number of cells in the source domain
ncP = 50; % number of cells in the patterning domain
diameter = 4.9; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 19.26; % mean gradient length [µm]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d;

mu_a = (diameter/2)^2*pi % mean cell area 

CV_a = [0.01:0.01:0.09 0.1:0.05:0.95 1:0.5:10]';

CV = 0.3;

nc = ncS + ncP;

lambda = NaN(length(CV), 1);

p = mu_p * ones(nc, 1);

p = random(logndist(mu_p, mu_p * CV), nc, 1);

a = random(logndist(mu_a, mu_a * CV), nc, 1);

l = 2 * sqrt(a / (pi));

x0 = (-ncS:ncP) * diameter;
x0 = sort([x0 x0(2:end-1)])

x = [1,2,3,4,5]

x(1:end-1)


% log-normal distribution with adjusted mean & stddev
function pd = logndist(mu, sigma)
    pd = makedist('Lognormal', 'mu', log(mu/sqrt(1+(sigma/mu)^2)), 'sigma', sqrt(log(1+(sigma/mu)^2)));
end
end
