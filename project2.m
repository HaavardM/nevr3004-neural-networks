close all; 
set(groot,'defaultLineLineWidth',2.0)
mkdir('report/project2/figs');

N = 50;
S = zeros(N,1);
noise = 0.2;
V = rand(N,1);
V(V >= 0.5) = 1;V(V < 0.5) = -1;
% Create matrix with w_ij = (v_i*v_j/N)
W = V * V';
% Reset diagonal elements to zero (no self connections)
assert(energy(V, W) == energy2(V, W));

% Task 1
T = 100;
[M, H, S, E] = runSim(V, W, V, T);
fig1 = figure(1);
plot(1:T, M); %The network is stable
xlabel('Iterations');
ylabel('Similarity');
title(["Stability of network", "with perfect initialization"]);
saveas(fig1, "report/project2/figs/stable.eps", "epsc");
% Task 2
V_noise = patternWithNoise(V, noise);
M = runSim(V_noise,W,V,T);
fig2 = figure(2);
plot(1:T, M); %The network is stable, and "reconstructs" the original pattern
xlabel('Iterations');
ylabel('Similarity');
title(["Stability of network", "with noisy initialization"]);

saveas(fig2, "report/project2/figs/stable-with-noise.eps", "epsc");

%% Task 3
U = rand(N, 1);
U(U >= 0.5) = 1; U(U < 0.5) = -1;
U_noise = patternWithNoise(U, noise);
W = (V*V' + U*U') / 2;
[M, H, S, E] = runSim([V,V_noise, U, U_noise], repmat(W, 1, 1, 4), [V, V, U, U], T);
fig3 = figure(3);
plot(1:T, M);
xlabel('Iterations');
ylabel('Similarity');
title('Multiple patterns');
l = legend('V vs V', 'V_{noise} vs V', 'U vs U', 'U_{noise} vs U');
l.Location = 'southeast';
saveas(fig3, "report/project2/figs/multiple-patterns.eps", "epsc");
fig4 = figure(4);
plot(1:T, E);
xlabel('Iterations');
ylabel('Energy');
title('Network energy');
legend("V", "V_{noise}", "U", "U_{noise}");
saveas(fig4, "report/project2/figs/multiple-patterns-energy.eps", "epsc");
%% QR-codes
noise = 0.9;
fig5 = figure(5);
I = loadQRPatternFromPng('qr-code.png');
V = I(:);
[I, N] = loadQRPatternFromPng('qr-code-2.png');
U = I(:);
V_noise = patternWithNoise(V, noise);
U_noise = U;
U_noise(floor(end/2):end) = 1;
U_inv = -1*U;
U_inv2 = U_inv;
U_inv2(1:floor(end/3)) = -1* U_inv2(1:floor(end/3));
R = rand(size(V,1), 1);
W = V*V' + U*U';
for i = 1:size(R, 2)
    W = W + R(:, i)*R(:, i)';
end
W = W / 4; %M=4
T = 1000;
[M, H, S, E] = runSim([V, V_noise, U, U_noise, U_inv, U_inv2], repmat(W, 1, 1, 6), [V, V, U, U, U, U], T);
subplot(2,3,1);
displayImage(V);
title('Original QR#1');
axis image;
subplot(2,3,2);
displayImage(V_noise);
title(sprintf('With %.0f %% noise', noise * 100));
axis image;
subplot(2, 3, 3);
displayImage(S(:, 2));
title(["Reconstructed", "by Hopfield network"]);
axis image;
subplot(2,3,4);
displayImage(U);
title('Original QR#2');
axis image;
subplot(2,3,5);
displayImage(U_noise);
title('After loosing half');
axis image;
subplot(2, 3, 6);
displayImage(S(:, 4));
title(["Reconstructed", "by Hopfield network"]);
axis image;
% Create arrow
annotation('arrow',[0.63 0.68],...
    [0.75 0.75]);

% Create arrow
annotation('arrow',[0.35 0.4],...
    [0.75 0.75]);

% Create arrow
annotation('arrow',[0.35 0.4],...
    [0.275 0.275]);

% Create arrow
annotation('arrow',[0.63 0.68],...
    [0.275 0.275]);
saveas(fig5, 'report/project2/figs/qr-code.eps', "epsc");
fig6 = figure(6);
plot(1:T, M(:, 1:4));
l = legend('QR#1', 'QR#1 with noise', 'QR#2', 'QR#2 with loss');
l.Location = 'southeast';
xlabel('Iterations');
ylabel('Similarity');
title("Stability of network");
saveas(fig6, 'report/project2/figs/qr-code-sim.eps', "epsc");

fig7 = figure(7);
plot(1:T, E(:, 1:4));
xlabel("Iterations");
ylabel("Energy");
title("Energy function for Hopfield network");
legend('QR#1', 'QR#1 with noise', 'QR#2', 'QR#2 with loss');
saveas(fig7, 'report/project2/figs/qr-code-energy.eps', "epsc");

fig8 = figure(8);
subplot(2,3,1);
displayImage(U);
title('Original QR#2');
axis image;
subplot(2,3,2);
displayImage(U_inv);
title('After invertion');
axis image;
subplot(2, 3, 3);
displayImage(S(:, 5));
title(["Reconstructed", "by Hopfield network"]);
axis image;
subplot(2,3,4);
displayImage(U);
title('Original QR#2');
axis image;
subplot(2,3,5);
displayImage(U_inv2);
title({'After inverting', 'a 2/3 of the bits'});
axis image;
subplot(2, 3, 6);
displayImage(S(:, 6));
title(["Reconstructed", "by Hopfield network"]);
axis image;
% Create arrow
annotation('arrow',[0.63 0.68],...
    [0.75 0.75]);

% Create arrow
annotation('arrow',[0.35 0.4],...
    [0.75 0.75]);

% Create arrow
annotation('arrow',[0.35 0.4],...
    [0.275 0.275]);

% Create arrow
annotation('arrow',[0.63 0.68],...
    [0.275 0.275]);
saveas(fig8, "report/project2/figs/qr-inverted.eps", "epsc");

function [V_noise] = patternWithNoise(V, p)
    V_noise = V;
    N = size(V, 1);
    n_p = floor(p*N);
    V_noise(randsample(N, n_p)) = rand(n_p, 1);
    V_noise(V_noise >= 0.5) = 1; V_noise(V_noise < 0.5) = -1;
end

function [I, N] = loadQRPatternFromPng(filename)
    I = rgb2gray(imread(filename));
    I = I(any(I == 0, 2), :);
    I = I(:, any(I == 0, 1));
    I = double(I);
    I(I < 50) = -1; 
    I(I > 50) = 1;
    N = numel(I);
end

function displayImage(I)
    N = sqrt(numel(I));
    img = reshape(I, N, N);
    img(img == 1) = 255;
    img(img == -1) = 0;
    img = uint8(img);
    image(img);
    yticks([]);
    xticks([]);
    colormap(gray);
end

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