function fit_line

% read the data 

% data where the cell diameter is varied and cv_k = 0.3 & cv_area = 0.5
T_diff = readtable('Parameter_tables/Mean_Diameter_change_CV_Diff.csv');
T_all = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_all.csv');
T_p = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_p.csv');
T_d = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_d.csv');


% data where the cell area variability is varied
T_diff_area = readtable('Parameter_tables/Cell_Area_CV_Diff_length_version.csv');
T_all_area = readtable('Parameter_tables/Cell_Area_vs_CV_length_version_all.csv');
T_p_area = readtable('Parameter_tables/Cell_Area_vs_CV_length_version_p.csv');
T_d_area = readtable('Parameter_tables/Cell_Area_vs_CV_length_version_p.csv');

% original publication where domain lenght L is changed, cv_area = 0.5
% added
T_length_all =  readtable('script_from_paper_cv_added/L_vs_CV_vs_cv_area_0.5_all.csv');

% Import data to find average cell diameter
T_area_data = readtable('Area_Data/average_cell_area_different_tissues.csv')

% Import data cv_area = 0.5 added 
T_cv_area_added_all = readtable('script_from_paper_cv_added/CV_vs_CV_cv_area_0.5_all.csv')
T_cv_area_added_p = readtable('script_from_paper_cv_added/CV_vs_CV_cv_area_0.5_p.csv')
T_cv_area_added_d = readtable('script_from_paper_cv_added/CV_vs_CV_cv_area_0.5_d.csv')
T_cv_area_added_D = readtable('script_from_paper_cv_added/CV_vs_CV_Diff_cv_area_0.5.csv')
%=========================================================================%
% model using polyfit()
%=========================================================================%

% get slope for CV_lambda 

CV_diff_slope = linearfit(log(T_diff.diameter_lambda), log(T_diff.CV_lambda))

m = CV_diff_slope(1)
k = CV_diff_slope(2)

CV_all_slope = linearfit(log(T_all.diameter_lambda), log(T_all.CV_lambda))

sub_p_cv_lambda = T_p.CV_lambda(1:10)

CV_p_slope = linearfit(log(T_p.diameter_lambda(1:10)), log(T_p.CV_lambda(1:10)))


% get slope for CV_c0

CV_diff_slope_c0 = linearfit(log(T_diff.diameter_lambda), log(T_diff.CV_0))
CV_all_slope_c0 = linearfit(log(T_all.diameter_lambda), log(T_all.CV_0))
CV_p_slope_c0 = linearfit(log(T_p.diameter_lambda), log(T_p.CV_0))
CV_d_slope_c0 = linearfit(log(T_d.diameter_lambda), log(T_d.CV_0))

function f = linearfit(x, y)
    f = polyfit(x,y,1);
end

%=========================================================================%
% fit nlm instead of polyfit 
%=========================================================================%

% set up the model 
beta0 = [0,0];
beta_q = [0,0,0];
linear_function = @(b,x)(b(1)*x + b(2));
quadratic_function = @(b,x)(b(1)*x.^2 + b(2)*x + b(3));

% model for diffusion coeff. 
mdl_diff = fitnlm(log(T_diff.diameter_lambda), log(T_diff.CV_lambda), linear_function, beta0)

%model for degradation d 
mdl_d = fitnlm(log(T_d.diameter_lambda), log(T_d.CV_lambda), linear_function, beta0)

% model for production p 
mdl_p = fitnlm(log(T_p.diameter_lambda), log(T_p.CV_lambda), linear_function, beta0)

% model with all parameters
mdl_all = fitnlm(log(T_all.diameter_lambda), log(T_all.CV_lambda), linear_function, beta0)

% fit the same models for CV_0
mdl_all_cv0 = fitnlm(log(T_all.diameter_lambda), log(T_all.CV_0), linear_function, beta0)
mdl_p_cv0 = fitnlm(log(T_p.diameter_lambda), log(T_p.CV_0), linear_function, beta0)
mdl_diff_cv0 = fitnlm(log(T_diff.diameter_lambda), log(T_diff.CV_0), linear_function, beta0)
mdl_d_cv0 = fitnlm(log(T_d.diameter_lambda), log(T_d.CV_0), linear_function, beta0)

% fit the domain lenght vary data with cv_area = 0.5 added 
mdl_length_cv_0 = fitnlm(T_length_all.L, T_length_all.CV_0, quadratic_function, beta_q)
mdl_length_cv_lambda = fitnlm(log(T_length_all.L), log(T_length_all.CV_lambda), linear_function, beta0)


%=========================================================================%
% Use area data to calculate average diameter of cells 
%=========================================================================%

names = T_area_data.(1)
areas = T_area_data.(2)

% calcualte diameters for each cell type
diameter = 2*sqrt(areas/pi)

%=========================================================================%
% Fit CV_area data
%=========================================================================%

CV_area = T_diff_area.CV_area
CV_area_const = CV_area(1:30);
CV_area_linear = CV_area(30:end);
area_diff = T_diff_area.CV_lambda;
constant_area_diff = area_diff(1:30);
linear_area_diff = area_diff(30:end);

mdl_const_area_diff = fitnlm(log(CV_area_const), constant_area_diff, linear_function, beta0)
mdl_const_area_diff = fitnlm(log(CV_area_linear), linear_area_diff, linear_function, beta0)

%=========================================================================%
% Fit CV_area = 0.5 added to initial script by Vetter & Iber 
%=========================================================================%

cv_lambda_all_initial_script_cv_area_added = fitnlm(log(T_cv_area_added_all.CV),log(T_cv_area_added_all.CV_lambda), quadratic_function, beta_q)

% CV_C0
subset_c0_p = T_cv_area_added_p.CV_0(1:22); 
subset_c0_d = T_cv_area_added_d.CV_0(1:22);
subset_c0_D = T_cv_area_added_D.CV_0(1:22);

subset_c0_d

cv_d = T_cv_area_added_all.CV(1:22);

subset_c0_all = T_cv_area_added_all.CV_0(1:22);

cv_0_all_initial_script_cv_area_added = fitnlm(cv_d,subset_c0_all, linear_function, beta0)
cv_0_p_initial_script_cv_area_added = fitnlm(cv_d,subset_c0_p, linear_function, beta0)
cv_0_d_initial_script_cv_area_added = fitnlm(cv_d,subset_c0_d, linear_function, beta0)
cv_0_D_initial_script_cv_area_added = fitnlm((cv_d),subset_c0_D, linear_function, beta0)

%=========================================================================%
% plot some data to check
%=========================================================================%

x= linspace(0.01,10, 1000);
plot(T_diff.diameter_lambda, T_diff.CV_lambda)
plot(x, x.^m.*exp(k))
set(gca, 'XScale', 'log', 'YScale', 'log')


end
