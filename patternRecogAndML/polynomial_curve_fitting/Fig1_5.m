#!/bin/octave

# training data
training_N = 10;
training_x = [0:training_N-1] / (training_N-1);
training_t = sin(2*pi*training_x) + .3*randn(1, training_N);

# test data
test_N = 100;
test_x = [0:test_N-1] / (test_N-1);
test_t = sin(2*pi*test_x) + .3*randn(1, test_N);

# functions
function ret = y(x, w)
	ret = x' .^ [0:length(w)-1] * w;
endfunction

function error = least_squares(x, t, w)
	sum = 0;
	N = length(x);
	for n = 1:N
		sum += (y(x(n), w) - t(n)) ^ 2;
	endfor
	error = sum;
endfunction


Max_M = 9;
RMS = zeros(Max_M+1, 2);
for M = 0:Max_M
	# Calculating coefficients vector w* of Mth degree polynomial so that dleast_squares/dw = 0 at w*
	sum_x_powers = zeros(1,2*M+1);
	for n=1:training_N
		sum_x_powers += training_x(n) .^ [0:2*M];
	endfor
	x_powers_matrix = zeros(M+1);
	for j=1:M+1
		x_powers_matrix([j*(M+1)-M:j*(M+1)]) += sum_x_powers(j:j+M);
	endfor
	target = zeros(M+1, 1);
	for n=1:training_N
		target += training_x(n) .^ [0:M]' * training_t(n);
	endfor

	w = inv(x_powers_matrix) * target;

	
	RMS(M+1, 1) = sqrt(2*least_squares(training_x, training_t, w)/training_N);
	RMS(M+1, 2) = sqrt(2*least_squares(test_x, test_t, w)/test_N);
	# printf("%d %d %d\n", M, RMS(M+1, 1), RMS(M+1, 2));
endfor

# plot
p = plot([0:9], RMS([1:Max_M+1]), "-ob", [0:9], RMS([Max_M+2:2*(Max_M+1)]), "-or");
axis([-1, 10, 0, 1.5]);
refresh()
waitfor(p)
