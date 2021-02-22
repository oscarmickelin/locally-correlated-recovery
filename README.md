The following files accompany the paper https://arxiv.org/abs/2102.09661

The code requires the Tensor Toolbox (https://www.tensortoolbox.org/) to be installed.

## Files:
### Example Files:

each of these run a different example

ex_bound_bs.m

ex_error_tol.m

ex_num_iter.m


### Main files:
solve_unknown_twosided.m - Sets up a randomly generated recovery instance and solves it using the algorithm of the paper

solve_unknown_noise_twosided.m - Sets up a noisy version of a randomly generated recovery instance and solves it using the algorithm of the paper

solve_first_twosided.m - Iteratively solves the linear system in step 2 of the algorithm

solve_second_twosided.m - Iteratively solves the linear system in step 3 of the algorithm

solve_third.m - Solves step 5 in the algorithm

solve_sub_slice_second_twosided.m - Updates one matrix X_num in step 3 of the algorithm

solve_sub_slice_twosided.m - Updates one matrix X_num in step 2 of the algorithm

solve_each_slice_twosided.m - solves one iteration of AX^T - XA^T, where X is a matrix with structure band-diagonal + cross

solve_each_slice_second_twosided.m - solves one iteration of AX - X^TA^T, where X is a matrix with band structure

jennrich.m - Implements Jennrich's algorithm for orthogonally decomposable tensors

calc_stop.m - Calculates the relative improvement between current iterate Xs and the previous iterate Xsprev, for elements contained in bad_rows
