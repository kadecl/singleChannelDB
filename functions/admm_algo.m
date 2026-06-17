function speech_restore_bin = admm_algo(reverb_speech_bin, rho, M, N)
Lambda_1 = zeros(N, M);
Lambda_2 = zeros(N, M);
lambda = zeros(N, 1);
Z = [reverb_speech_bin zeros(N, M-1)];
Y_1 = zeros(N, M);
C_circ = zeros(N, M);
R = zeros(N, M);
Z_circ = zeros(N, M);
l = 0;

while l < 1000
    l = l + 1;

    Y_1_hat = Z - Lambda_1;
    for n = 1 : N
        Y_1(n,:) = max(1 - 1 ./ (rho * norm(Y_1_hat(n,:))),0) * Y_1_hat(n,:);
    end

    Q_hat = Z - Lambda_2;
    [U_1,S_1,V_1] = svds(Q_hat, 1);
    Y_2 = U_1 * S_1 * V_1';

    C = Y_1 + Lambda_1 + Y_2 + Lambda_2;
    for c = 1 : M
        C_circ(:, c) = circshift([C(:,c)],c-1);
    end
    C_circ_1 = sum(C_circ, 2);
    w = 2 .* (reverb_speech_bin + lambda) - C_circ_1;
    for c = 1 : M
        R(:, c) = circshift(w./(M+2), -c + 1);
    end
    Z = (C + R) ./ 2;
    for c = 1 : M
        Z_circ(:, c) = circshift([Z(:,c)], c - 1);
    end
    Z_circ = sum(Z_circ, 2);
    Lambda_1 = Lambda_1 - Z + Y_1;
    Lambda_2 = Lambda_2 - Z + Y_2;
    lambda = lambda - Z_circ + reverb_speech_bin;
end

[U_2,S_2,V_2] = svds(Z, 1);
ir_restore_bin = V_2(:,1);
speech_restore_bin = S_2(1, 1) * U_2(:, 1) ./ sign(ir_restore_bin(1));

end