clear all;
close all;
clc;

x = linspace(0, 1, 11);

data.eps1 = 0.5;
data.eps2 = 1;

i = 1;

if (i == 1)
    # Non-homogeneous.
    c = -0.4375 * data.eps2 / (0.5 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    d = 1.75 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) / (sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + 2 * data.eps2 * sinh(0.5 / sqrt(data.eps1)));

    u1 = 1 + 2 * c * sinh(x / sqrt(data.eps1));
    u2 = -0.5 * (x - 1) .* (x + 2 * d);
else
    # Homogeneous.
    c = -0.5 * data.eps2 / (0.5 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    d = 2 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) / (sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + 2 * data.eps2 * sinh(0.5 / sqrt(data.eps1)));

    u1 = 1 + 2 * c * sinh(x / sqrt(data.eps1));
    u2 = -d * (x - 1);
end

plot(x, u1, x, u2);
legend('u_1', 'u_2');
