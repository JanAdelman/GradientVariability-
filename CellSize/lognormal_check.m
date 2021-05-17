T = readtable('figure1_coordinates_ventral.csv');

% fit log normal to ventral data of figure 1
parmhat = lognfit(T.x);

parmhat

% try with second approach 
pd = fitdist(T.x,'Lognormal')

pd

hold on

% plot the data 
scatter(T.x, T.y)

x_values = 0:1:60;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
plot(pd,y)
h = gca;

% it does not fit the distribution