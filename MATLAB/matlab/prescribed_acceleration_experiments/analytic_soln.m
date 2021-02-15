function f = analytic_soln(x, t)

string_parameters

n_terms = 20;
disp_omega = 2*pi*330;

a = zeros(n_terms,1);
b = zeros(n_terms,1);
for n = 1:2:n_terms
    kn = n*pi/len;
    a(n) = -4*disp_omega^2/(n*pi*(disp_omega^2 - (kn*c)^2));
    b(n) =  4*disp_omega*kn*c/(n*pi*(disp_omega^2 - (kn*c)^2));
end

f = ones(size(x)) * sin(disp_omega*t);

for n=1:length(a)
    kn = n*pi/len;
    f = f + a(n)*sin(kn*x)*sin(disp_omega*t);
    f = f + b(n)*sin(kn*x)*sin(kn*c*t);
end
