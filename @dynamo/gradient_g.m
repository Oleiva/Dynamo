function ret = gradient_g(self, t, k, c)
% Gradient of the auxiliary function g(A,B), at
% time slice t, ensemble member k, with respect to u_c.
%
% g(A,B) = \tr_S(A'*B) is enough to compute the error function
% E(A,B) whenever |B| is constant. This holds in every closed
% system task except "closed state_partial".
%
% Actually returns the gradient times -V*U', for the benefit of E_abs.

% Uses H{t} and L{t+1}, and additionally
% U{t_c}, U{t_tau+1} (fd, eig, auxmatrix) or
% U{t+1} (series).


if c < 0
    % tau control:  dP/dtau = G P = P G  (exact)
    ret = self.seq.tau_deriv(t) * trace_matmul(self.cache.L{t+1, k}, self.cache.H{t, k} * self.cache.U{t+1, k});
else
    % other control
    switch self.config.dP
      case 'fd'
        % f'(x) = (f(x + eps) -f(x))/eps
        % Trivial and relatively slow, but a good reference point.
        [P, epsilon] = self.finite_diff_P(t, k, c);
        % We may compute the finite diff approximation at three locations: P, g, E
        %dPdu = (P -self.cache.P{t, k}) / epsilon;
        %ret = -self.cache.VUdagger * trace_matmul(self.cache.L{t+1, k}, dPdu * self.cache.U{t, k});
        %return
        g = trace_matmul(self.cache.L{t+1, k}, P * self.cache.U{t, k});
        %ret = -(self.cache.VUdagger / epsilon) * (g -self.cache.g{k});
        %return
        E = self.config.f_max -abs(g);
        %E = self.config.f_max -real(g);
        ret = (E -self.cache.E) / epsilon;
        return

      case 'eig'
        dPdu = dPdu_eig(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});
        ut = t;

      case 'aux'
        dPdu = dPdu_auxmatrix(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c});
        ut = t;

      case 'series'
        % This test is _really_ expensive, around 20% of total running time.
        %self.gradient_test(t, k, c);
        dPdu = dPdu_series(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c}, 10);
        %dPdu = self.system.B{k, c};  % series, 1st order
        ut = t+1;
    end
    ret = self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.L{t+1, k}, dPdu * self.cache.U{ut, k});
end

ret = -ret * self.cache.VUdagger;
end
