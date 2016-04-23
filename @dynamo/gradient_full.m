function ret = gradient_full(self, t, k, c)
% Gradient of E_full(A,B), at time slice t, ensemble member k, with respect to u_c.

% Uses H{t} and L{t+1}, and additionally
% U{t_c}, U{t_tau+1} (eig, auxmatrix, fd, series_ss) or
% U{t+1} (series)

% \tr_E(B) -A
delta = self.cache.g{k} -self.system.X_final;

if c < 0
    % tau control:  dP/dtau = G P = P G  (exact)
    ret = self.seq.tau_deriv(t) * inprod(delta, partial_trace(self.cache.L{t+1, k} * self.cache.H{t, k} * self.cache.U{t+1, k}, self.system.dimSE, 2));
else
    % other control
    switch self.config.dP
      case 'fd'
        % f'(x) = (f(x + eps) -f(x))/eps
        % Trivial and relatively slow, but a good reference point.
        [P, epsilon] = self.finite_diff_P(t, k, c);
        % We may compute the finite diff approximation at three locations: P, B_S, E
        %dPdu = (P -self.cache.P{t, k}) / epsilon;
        %dB_Sdu = partial_trace(self.cache.L{t+1, k} * (dPdu * self.cache.U{t, k}), self.system.dimSE, 2);
        %ret = inprod(delta, dB_Sdu);
        %return
        B_S = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 2);
        %dB_Sdu = (B_S -self.cache.g{k}) / epsilon;
        %ret = inprod(delta, dB_Sdu);
        %return
        E = 0.5 * norm2(B_S -self.system.X_final);
        ret = (E -self.cache.E) / epsilon;
        return

      case 'eig'
        dPdu = dPdu_eig(self.cache.H_v{t, k}, self.cache.H_eig_factor{t, k}, self.system.B{k, c});
        ut = t;

      case 'aux'
        dPdu = dPdu_auxmatrix(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c});
        ut = t;

      case 'series_ss'
        %self.gradient_test(t, k, c);
        W = self.cache.W{t, k};
        n = length(W);  % number of squarings used for P
        s = 2^n;  % corresponding scale
        G = (self.seq.tau(t) / s) * self.cache.H{t, k};  % scaling
        % temp == expm(G)
        if n == 0
            temp = self.cache.P{t, k};
        else
            temp = W{1};
        end
        % d(expm(G))/du,  NOTE missing tau here
        dPdu = dPdu_series(G, self.system.B{k, c} / s, 7) * temp;
        % squaring
        for r=1:n
            % W{r} == expm(G)^(2^(r-1))
            dPdu = dPdu * W{r} +W{r} * dPdu;  % chain rule for the derivative
        end
        ut = t;  % no P{t} required any more

      case 'series'
        % NOTE missing tau here
        dPdu = dPdu_series(self.seq.tau(t) * self.cache.H{t, k}, self.system.B{k, c}, 12);
        ut = t+1;  % provide the missing P{t}
    end
    % NOTE adding the missing tau factor
    ret = self.seq.tau(t) * self.seq.fields_deriv(t, c) * inprod(delta, partial_trace(self.cache.L{t+1, k} * dPdu * self.cache.U{ut, k}, self.system.dimSE, 2));
end
end
