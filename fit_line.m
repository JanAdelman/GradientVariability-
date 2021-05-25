function fit_line

T_diff = readtable('Parameter_tables/Mean_Diameter_change_CV_Diff.csv');
T_all = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_all.csv');
T_p = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_p.csv');
T_d = readtable('Parameter_tables/Mean_Diameter_change_vs_CV_d.csv');


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

x= linspace(0.01,10, 1000);
plot(T_diff.diameter_lambda, T_diff.CV_lambda)
plot(x, x.^m.*exp(k))
set(gca, 'XScale', 'log', 'YScale', 'log')


end
