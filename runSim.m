function [M, H, S, E] = runSim(S, W, V, T)
    P = size(V, 2);
    N = size(S, 1);
    M = zeros(T, P);
    H = zeros(N, P);
    E = zeros(T, P);
    
    for p = 1:P
        it = randperm(N);
        for t = 1:T
            W_t = W(:, :, p);
            W_t(1:N+1:end) = 0;
            E(t, p) = energy2(S(:, p), W_t);
            M(t, p) = S(:, p)' * V(:,p)  / N;
            H(:, p) = W_t * S(:, p);
            n = it(mod(t, N) + 1);
            S(n, p) = sign(H(n, p));
        end
    end
end

function E = energy(S, W)
    x = ones(numel(S), 1);
    E = (-1/2) * x'* (W.*(S*S'))*x;
end

function E = energy2(S, W)
    E = (-1/2)*S'*W*S;
    assert(E == energy(S, W));
end