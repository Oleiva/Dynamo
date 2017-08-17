function term_reason = search(self, varargin)
% Run the optimization.


%% MATLAB-style options processing.
% Converts a list of fieldname, value pairs in varargin to a struct.
user_options = struct(varargin{:});

if isempty(self.opt)
    % no previous optimization rounds, this is the first one
    % set default termination conditions and other options
    self.opt.options = struct(...
        'control_mask',      self.seq.mask,...  % default mask
        'error_goal',        0.5 * (1e-4)^2 / self.system.norm2,...
        'max_evals',         1e6,...
        'max_walltime',      1800,...
        'max_cputime',       1e6,...
        'min_gradient_norm', 1e-20,...
        'plot_interval',     1);   % how often should we plot intermediate results?
    self.opt.matlab_options = struct();
else
    % there has been at least one optimization round before
    % copy all the options from the last round, modify them with the user input.
end

% initialization of the optimization data structures
% self.opt.options and self.opt.matlab_options are kept as is.
self.init_opt();

% modify self.opt.options with user_options (if any)
[self.opt.options, unused] = apply_options(self.opt.options, user_options, true);
% the rest are dumped into matlab_options
[self.opt.matlab_options] = apply_options(self.opt.matlab_options, unused, false);


%% run the optimizer

fprintf('Optimization space dimension: %d\n', sum(sum(self.opt.options.control_mask)));

% define the optimization problem
obj_func = @(x) goal_and_gradient_function_wrapper(self, x);

% run BFGS optimization
[self.opt.matlab_exitflag, self.opt.matlab_output] = self.search_BFGS(obj_func, self.opt.matlab_options);


%% make a copy of the current run's statistics

self.stats{end+1} = self.opt;


%% termination reason

if self.opt.matlab_exitflag == -1
    % term_reason comes from monitor_func
else
    % terminated by optimizer
    self.opt.term_reason = self.opt.matlab_output.message;
end
term_reason = self.opt.term_reason;
end


function [err, grad] = goal_and_gradient_function_wrapper(self, x)
% x is a vector containing (a subset of) the controls

    self.opt.n_eval = self.opt.n_eval +1;

    self.update_controls(x, self.opt.options.control_mask);
    [err, grad] = self.compute_error(self.opt.options.control_mask);
    self.opt.last_grad_norm = sqrt(sum(sum(grad .* grad)));
end
