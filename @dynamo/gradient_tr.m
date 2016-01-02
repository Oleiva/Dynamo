function ret = gradient_tr(self, t, k, c)
% Gradient of error_tr.

% Uses H{t} and L{t+1}, and additionally
% U{t_c}, U{t_tau+1} (fd, eig, auxmatrix) or
% U{t+1} (series).

if c < 0
    % dP_t/dtau_{t} = H_t P_t = P_t H_t
    dQdtau = partial_trace(self.cache.L{t+1, k} * self.cache.H{t, k} * self.cache.U{t+1, k}, self.system.dimSE, 1);
    ret = -self.seq.tau_deriv(t) * trace_matmul(self.cache.VUdagger, dQdtau);
else
    % other control
    switch self.config.dP
      case 'fd'
        % f'(x) = (f(x + eps) -f(x))/eps
        % Trivial and relatively slow, but a good reference point.
        [P, epsilon] = self.finite_diff_P(t, k, c);
        % We may compute the finite diff approximation at three locations: P, g, E
        %dPdu = (P -self.cache.P{t, k}) / epsilon;
        %dgdu = partial_trace(self.cache.L{t+1, k} * (dPdu * self.cache.U{t, k}), self.system.dimSE, 1);
        g = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 1);
        %dgdu = (g -self.cache.g{k}) / epsilon;
        %ret = -trace_matmul(self.cache.VUdagger, dgdu);
        %return
        E = self.config.f_max -sum(svd(g));
        ret = (E -self.cache.E) / epsilon;

      case 'eig'
        dPdu = dPdu_eig(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});
        dQdu = partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{t, k}, self.system.dimSE, 1);
        ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.VUdagger, dQdu);
    end
end
