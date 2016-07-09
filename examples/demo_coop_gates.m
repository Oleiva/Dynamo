function dyn = demo_coop_gates()
% Example: Concurrent optimization method of co-operative gates.
%
% Assume we wish to perform a Ramsey sequence:
% Start in the |0> state, do a pi/2_y rotation to the Bloch sphere
% equator, wait for a time t while the state gains a relative phase,
% do a -pi/2_y rotation, measure |0><0|.
%
% For the sequence to work, the rotations do not have to be exact y rotations.
% It's enough if they "co-operate" by bringing the state to the
% equator and back with phases that cancel each other. This extra
% degree of freedom can make gates shorter and/or easier to optimize.
%
%! M. Braun and S. Glaser, New J. Phys. 16, 115002 (2014).

% Ville Bergholm 2015-2016


% Pauli Z
Z = diag([1, -1]);

% Liouvillian superoperator, projects to the equator of the Bloch sphere.
P_xy = [1, 0, 0, 1; 0, 2, 0, 0; 0, 0, 2, 0; 1, 0, 0, 1] / 2;

dim = 2;  % one qubit


%% Set up Dynamo

% Idea: optimize a state transfer sequence from |0> to |0>, but in
% the middle of the sequence project the state to the equator.
% Unless the first half of the sequence has already brought the
% state to the equator, it will lose some purity and thus the
% optimization can no longer reach zero error.
% The result of this optimization is a sequence where the first and
% second halves represent the two co-operative pi/2 rotations.

task = 'open state overlap';
ini = [1, 0].';
fin = ini;
H_drift = Z / 2;
[H_ctrl, c_labels] = control_ops(dim, 'xy');

T = 2;
n_bins = 7;

dyn = dynamo(task, ini, fin, H_drift, H_ctrl);
dyn.system.set_labels('Co-operative gates optimization demo', dim, c_labels);
dyn.seq_init(n_bins, T * [0.5, 1.5]); %, control_type, control_par);

% random, constant initial controls
dyn.set_controls(0.1 * randn(1, 2));

mask = dyn.full_mask(false);

t = ceil(n_bins/2);  % replace the propagator in this bin by P_xy
dyn.cache.set_P(t, 1, P_xy);
mask(t,:) = false;  % do not update the controls so as not to overwrite P_xy
dyn.seq.raw(t, 1:end-1) = 0; % set the phantom controls to zero to make the sequence look nicer

dyn.ui_open();
dyn.search(mask);
