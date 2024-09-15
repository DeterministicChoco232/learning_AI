#!/bin/octave

# training data
N = 10;
x = [0:N-1] / (N-1);
t = sin(2*pi*x) + .3*randn(1, N);

# functions
function ret = y(x, w)
	ret = x' .^ [0:length(w)-1] * w;
endfunction

function error = modified_least_squares(x, t, w, lambda)
	sum = w' * w * lambda / 2;
	N = length(x);
	for n = 1:N
		sum += (y(x(n), w) - t(n)) ^ 2;
	endfor
	error = sum;
endfunction


for lambda = [exp(-18) exp(0)]
	M=9
	# Calculating coefficients vector w* of Mth degree polynomial so that dmodified_least_squares/dw = 0 at w*
	sum_x_powers = zeros(1,2*M+1);
	for n=1:N
		sum_x_powers += x(n) .^ [0:2*M];
	endfor
	x_powers_matrix = zeros(M+1);
	for j=1:M+1
		x_powers_matrix([j*(M+1)-M:j*(M+1)]) += sum_x_powers(j:j+M);
	endfor
	x_powers_matrix += lambda * eye(M+1);
	target = zeros(M+1, 1);
	for n=1:N
		target += x(n) .^ [0:M]' * t(n);
	endfor

	w = inv(x_powers_matrix) * target;



	printf("Modified least squares: %d\n", modified_least_squares(x, t, w, lambda));
	# plot
	p_x = [0:100] / 100;
	p = plot(p_x, sin(2*pi*p_x), "g", x, t, "ob", p_x, y(p_x, w), "r");
	axis([-0.04, 1.04, -2, 2]);
	refresh()
	waitfor(p)
endfor
