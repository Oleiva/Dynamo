function analyze(self)
% Analyzes the results of an optimization run.


err = self.compute_error();
fprintf('Final normalized error: %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	err, self.opt.wall_time(end), self.opt.cpu_time(end), self.opt.term_reason);

fprintf('Number of gradient evaluations: %d\n', self.opt.n_eval);
fprintf('Final sequence duration: %g\n', sum(self.seq.tau));


% rough error estimates
%d = real(eig(self.system.A));
%n = length(d);
%T = sum(self.seq.tau);
%e_max = 1-exp(-T*sum(d)/n)
%e_min = 1-sum(exp(-T*d))/n

% plot the final sequence and some analytics
figure()
ax = subplot(3, 1, 1);
self.plot_seq(ax);


%% plot the error

ax = subplot(3, 1, 2);
set_plotstyle(ax);
offset = 0;
for k=1:length(self.stats)
    semilogy(ax, self.stats{k}.wall_time+offset, abs(self.stats{k}.error));
    hold on;
    semilogy(ax, self.stats{k}.wall_time(end)+offset, abs(self.stats{k}.error(end)), 'ko');
    offset = offset +self.stats{k}.wall_time(end);
end
grid on;
title('Optimization error')
xlabel('wall time (s)')
ylabel(ax, 'normalized error')


%% plot control integrals

ax = subplot(3, 1, 3);
set_plotstyle(ax);
%set(ax, 'LineStyleOrder','--')
offset = 0;
for k=1:length(self.stats)
    plot(ax, self.stats{k}.wall_time+offset, self.stats{k}.control_integral);
    offset = offset +self.stats{k}.wall_time(end);
    hold on;
end
grid on;
title('Control integral')
xlabel('wall time (s)')
ylabel(ax, 'control integral')
end
