close all; 
set(groot,'defaultLineLineWidth',2.0)
mkdir('report/project2/figs');


N = 500;
S = zeros(N,1);
noise = 0.7;
V = rand(N,1);
V(V >= 0.5) = 1;V(V < 0.5) = -1;
% Create matrix with w_ij = (v_i*v_j/N)
W = V * V' / N;
% Reset diagonal elements to zero (no self connections)


% Task 1
T = 1000;
M = runSim(V, W, V, T);
figure(1);
plot(1:T, M); %The network is stable

% Task 2
V_noise = patternWithNoise(V, noise);
M = runSim(V_noise,W,V,T);
figure(2);
plot(1:T, M); %The network is stable, and "reconstructs" the original pattern

% Task 3
U = rand(N, 1);
U(U >= 0.5) = 1; U(U < 0.5) = -1;
U_noise = patternWithNoise(U, noise);
W = V*V' + U*U';
M = runSim([V,V_noise, U, U_noise], repmat(W, 1, 1, 4), repmat(U, 1, 4), T);
figure(3);
plot(1:T, M);
legend('V', 'V_{noise}', 'U', 'U_{noise}');
%The overlap deivates away from zero

I = loadQRFromPng('qr-code.png');
V = I(:);
V_noise = patternWithNoise(V, noise);
W = V*V';
T = 4000;
[M, H, S] = runSim([V, V_noise], repmat(W, 1, 1, 2), [V, V], T);
figure(4);
subplot(1,3,1);
displayImage(V);
title('Original');
axis image;
subplot(1,3,2);
displayImage(V_noise);
title(sprintf('With %.0f %% noise', noise * 100));
axis image;
subplot(1, 3, 3);
displayImage(S(:, 2));
title(["Reconstructed", "by Hopfield network"]);
axis image;
print('report/project2/figs/qr-code', '-depsc');
figure(5);
plot(1:T, M);
l = legend('V', 'V_{noise}');
l.Location = 'southeast';
print('report/project2/figs/qr-code-sim', '-depsc');



function [V_noise] = patternWithNoise(V, p)
    V_noise = V;
    N = size(V, 1);
    n_p = floor(p*N);
    V_noise(randsample(N, n_p)) = rand(n_p, 1);
    V_noise(V_noise >= 0.5) = 1; V_noise(V_noise < 0.5) = -1;
end

function I = loadQRFromPng(filename)
    I = rgb2gray(imread(filename));
    I = I(any(I == 0, 2), :);
    I = I(:, any(I == 0, 1));
    I = double(I);
    I(I < 50) = -1; 
    I(I > 50) = 1;
    
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

function [M, H, S] = runSim(S, W, V, T)
    P = size(V, 2);
    N = size(S, 1);
    M = zeros(T, P);
    H = zeros(N, P);
    for p = 1:P
        for t = 1:T
            W_t = W(:, :, p);
            W_t(1:N+1:end) = 0;
            
            M(t, p) = S(:, p)' * V(:,p)  / N;
            H(:, p) = W_t * S(:, p);
            n = randi(N, 1);
            S(n, p) = sign(H(n, p));
        end
    end
end