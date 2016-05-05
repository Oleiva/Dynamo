function dyn = demo_ensemble(n_ensemble, delta)
% Example: Ensemble optimization.
%
%  dyn = demo_ensemble(n_ensemble, delta)
%
%  Optimizes a robust single-qubit control sequence for the state transfer |0> \to |1>
%  over an ensemble of systems with slightly different detunings.
%  Uses RWA, ignores the fast-rotating term.

% Ville Bergholm 2014-2016

if nargin < 2
    delta = 1;
    if nargin < 1
        n_ensemble = 5;
    end
end

%% Pauli matrices etc.

Z = diag([1, -1]);


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
H_drift = -delta * Z / 2;

% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'xy');

% limit driving Rabi frequency to the interval [-1, 1]
control_type = 'mm';
pp = [-1, 2];
control_par = {pp, pp};

% pure state transfer
initial = [1 0]';
final = [0 1]';

% define the ensemble
detuning = linspace(-1, 1, n_ensemble);
% gaussian distribution of weights
weight = exp(-(detuning * 2).^2);
weight = weight / sum(weight);

dyn = dynamo('closed ket', initial, final, @(k) detuning(k) * H_drift, H_ctrl, weight);
dyn.system.set_labels('Single-qubit ensemble optimization demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = pi * 7/3;  % enough for the short CORPSE sequence
dyn.seq_init(201, T * [1, 0], control_type, control_par);
dyn.set_controls(randn(1, 2));


%% Now do the actual search

dyn.ui_open();
dyn.search();
%dyn.analyze();


%% Plot the sequence error over different detunings
% TODO make evaluating the error of a sequence over various systems easier/more efficient
% TODO compute error/gradient for just a single ensemble member

det = linspace(-2,2,100);
err = [];
for k=1:length(det)
    temp = dynamo('closed ket', initial, final, det(k) * H_drift, H_ctrl);
    temp.seq = dyn.seq;
    temp.cache_init();
    err(k) = temp.compute_error();
end
figure
xlabel('detuning');
ylabel('error');
grid on;
hold on;
plot(det, err);
plot(detuning, min(err), 'ko')
legend('sequence error', 'ensemble')
end
