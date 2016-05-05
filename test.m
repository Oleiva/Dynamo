function test(d)
% Unit tests for DYNAMO.
%
%  d is a Dynamo instance containing the optimization problem used.
%  If no d is given, uses one of the test suite problems.
    
% Ville Bergholm 2015


% tolerance for numerical errors
tol = 1e-10;

if nargin < 1
    d = test_suite(1);
end

n_timeslots = d.seq.n_timeslots();
X = d.X();
temp = d.seq.tau;


%% export, then import a sequence

in_polar = false;

% no TU defined, do not try to use it
[A, desc] = d.export(in_polar, false);
desc
d.import(A, in_polar);
assert_equal(d.seq.tau, temp, tol);
assert_equal(d.X(), X, tol);

d.system.set_TU(1e-3);
% TU defined, get rid of it
[A, desc] = d.export(in_polar, false);
desc
d.import(A, in_polar);
assert_equal(d.seq.tau, temp, tol);
assert_equal(d.X(), X, tol);

% TU defined, keep it
[A, desc] = d.export(in_polar, true);
desc
d.import(A, in_polar, 1e-3);
assert_equal(d.seq.tau, temp, tol);
assert_equal(d.X(), X, tol);



%% split some random bins, this must not change the effect of the sequence

bins = find(rand(1, n_timeslots) < 0.5);
n = randi([2, 5]);
d.split(bins, n);
assert_equal(d.X(), X, tol);



disp('All tests passed.');
end


function assert_equal(a, b, tol)
%  Are a and b equal up to tolerance tol?

assert(norm(a-b) <= tol);
end
