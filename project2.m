close all; 
set(groot,'defaultLineLineWidth',3.0)
mkdir('report/project2/figs');

col = autumn(5);
set(gca, 'colororder', col);

N = 50;
S = zeros(N,1);
noise = 0.2;
V = rand(N,1);
V(V >= 0.5) = 1;V(V < 0.5) = -1;
% Create matrix with w_ij = (v_i*v_j/N)
W = V * V';

% Reset diagonal elements to zero (no self connections)

% Task 1
T = 50;
fig1 = figure(1);
V_noise = patternWithNoise(V, noise);
[M, H, S, E] = runSim([V, V_noise],repmat(W, 1, 1, 2),[V, V],T);
plot(1:T, M*100); %The network is stable, and "reconstructs" the original pattern
xlabel('Iterations');
ylabel('Similarity [%]');
title(["Stability of network", "with single pattern", "compared to V"]);
l = legend('V', 'V_{noise}');
l.Location = 'southeast';
ylim([-110 110]);

saveas(fig1, "report/project2/figs/stable.eps", "epsc");
fig2 = figure(2);
plot(1:T, E);
legend('V', 'V_{noise}');
xlabel('Iterations');
ylabel('Energy');
title(["Energy of network", "with single pattern"]);
saveas(fig2, "report/project2/figs/stable-energy.eps", "epsc");

%% Task 3 - Multiple patterns
U = rand(N, 1);
U(U >= 0.5) = 1; U(U < 0.5) = -1;
U_noise = patternWithNoise(U, noise);
W = (V*V' + U*U') / 2;
[M, H, S, E] = runSim([V,V_noise, U, U_noise], repmat(W, 1, 1, 4), [V, V, V, V], T);
fig3 = figure(3);
plot(1:T, M*100);
xlabel('Iterations');
ylabel('Similarity [%]');
title({'Similarity', 'with multiple patterns', 'compared to pattern V'});
l = legend('V', 'V_{noise}', 'U', 'U_{noise}');
l.Location = 'southeast';
ylim([-110 110]);
saveas(fig3, "report/project2/figs/multiple-patterns.eps", "epsc");
fig4 = figure(4);
plot(1:T, E);
xlabel('Iterations');
ylabel('Energy');
title({'Energy of network', 'with multiple patterns'});
legend("V", "V_{noise}", "U", "U_{noise}");
saveas(fig4, "report/project2/figs/multiple-patterns-energy.eps", "epsc");

%% Task 3 - capacity limits
N = 50;
MEM = 10;
Cit = 100;
noise = linspace(0.1, 0.5, 5);
succ = zeros(MEM, length(noise));
for it = 1:Cit
for m = 1:MEM
for n = 1:length(noise)
U = rand(N,m);
U(U >= 0.5) = 1; U (U < 0.5) = -1;
U_noise = patternWithNoise(U, noise(n));
T = 2*N;
W = zeros(N, N);
for i = 1:m 
    W = W + U(:, i)*U(:, i)';
end
W = W / m - eye(N);
[M, H, S, E] = runSim(U_noise, repmat(W, 1, 1, m), U, T);
succ(m, n) = succ(m, n) + mean(M(end, :) == 1);
end
end
end
fig31 = figure(31);
plot(1:MEM, succ / Cit)
xlabel("Number of stored memories");
ylabel("Success rate [%]");
title(["Memory restore success rate", "for different noise levels"])
legend(compose("Noise: %.1f %%", noise*100));
ylim([0, 1.1]);
saveas(fig31, 'report/project2/figs/capacity.eps', "epsc");

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
W = W / 2; %M=4
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
plot(1:T, M(:, 1:4)*100);
l = legend('QR#1', 'QR#1 with noise', 'QR#2', 'QR#2 with loss');
l.Location = 'southeast';
xlabel('Iterations');
ylabel('Similarity [%]');
ylim([-110 110]);
title(["Similarity of QR-codes", "compared with original"]);
saveas(fig6, 'report/project2/figs/qr-code-sim.eps', "epsc");

fig7 = figure(7);
plot(1:T, E(:, 1:4));
xlabel("Iterations");
ylabel("Energy");
title(["Energy of network", "with multiple QR-codes"]);
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
    
    for i = 1:size(V, 2)
        V_noise(randsample(N, n_p), i) = rand(n_p, 1);
    end
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

function E = lowerEnergy(N) 
    E = -(1/2)*N*(N-1);
end