function dPdu = dPdu_series(Gtau, B, n)
% Compute the derivative dP/du_c using the commutator series method.
% For efficiency this function actually returns (dP/du_c) * P^{-1} / tau.

    dPdu = B;
    for k=2:n
        % commutator
        B = (Gtau*B -B*Gtau)/k;
        dPdu = dPdu +B;
    end
end
