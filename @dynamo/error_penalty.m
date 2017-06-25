function ret = error_penalty(self, k, t, c)
% Forbidden state penalty when using error_full.
%
%  err  = error_penalty(self, k)
%    Total penalty in ensemble member k.
%
%  grad = error_penalty(self, k, t, c)
%    Gradient of total penalty in ensemble member k, with respect to
%    control field c, at time slice t.


if nargin == 2
    ret = 0;
    % compute the penalty after every time slice
    for t=1:self.seq.n_timeslots()+1
        ret = ret +self.system.penalty * self.cache.U{t, k};
    end
    ret = real(ret);
    self.cache.pen = ret;
else
    % compute the gradient of the forbidden states penalty
    % TODO gradient_setup needs an update for the P:s used here!

    temp = self.cache.dPdu_tail;
    ret = self.system.penalty * temp;
    for w = t+1:self.seq.n_timeslots()
        temp = self.cache.P{w,k} * temp;
        ret = ret +self.system.penalty * temp;
    end
    % multiply by the proper scale factor given by the partial
    % derivative of the control transformation
    ret = self.cache.dPdu_scale * ret;
    % get rid of small numerical errors
    ret = real(ret);
end
