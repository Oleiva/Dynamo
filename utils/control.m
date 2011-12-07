function [C, control_type] = control(dim, ctrl, s)
% Returns a cell vector of angular momentum control operators.
%  C = control(dim, ctrl [, s])
%
%  dim is the system dimension vector.
%  ctrl is a string consisting of chars 'x', 'y' and 'z', denoting
%  which controls are to be generated (for all the subsystems).
%  Control operators are generated for the first s subsystems.

% Ville Bergholm 2011


n = length(dim);
if nargin < 3
    s = n;
end
  
ctrl = unique(lower(ctrl));
c = length(ctrl);
C = cell(s, c);
for k=1:s
    J = angular_momentum(dim(k));
    for j=1:c
        temp = ctrl(j) - 'x' + 1; % MATLAB indexing
        C{k,j} = mkron(speye(prod(dim(1:k-1))), J{temp}, speye(prod(dim(k+1:end))));
    end
end
C = C(:);

if nargout == 2
  control_type = char('.' + zeros(size(C)));
end
