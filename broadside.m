clear all; close all; clc

x = audioread('broadside.wav')'; 

fftsize = 512
mic_dist = 0.11
sound_velocity = 343
fs = 16000
load win.mat; % window used in STFT
gamma_noise = sinc(2*[1:fftsize/2-1]/fftsize*fs*mic_dist/sound_velocity); % theoretical coherences of diffuse noise, will be used in enhanced PHAT
delay_search_range = 10
tau_true = 0

% estimated noise statistics, will be used in GEVD-MUSIC
rn1 = zeros(1, fftsize/2-1);
rn2 = zeros(1, fftsize/2-1);
cn12 = zeros(1, fftsize/2-1);
t = 1; fn = 0;
while t+fftsize-1 < 199957
    X = fft(win.*x(:, t:t+fftsize-1), [], 2);
    rn1 = rn1 + X(1, 2:fftsize/2).*conj(X(1, 2:fftsize/2)); % the first bin is DC, thus discarded
    rn2 = rn2 + X(2, 2:fftsize/2).*conj(X(2, 2:fftsize/2));
    cn12 = cn12 + X(1, 2:fftsize/2).*conj(X(2, 2:fftsize/2));
    t = t + 160;
    fn = fn + 1;
end
rn1 = rn1/fn; rn2 = rn2/fn; cn12 = cn12/fn;
gamma_noise_hat = (cn12)./sqrt((rn1).*(rn2));
figure;
plot([1:fftsize/2-1]/fftsize*fs, real(gamma_noise_hat), 'b-')
hold on; plot([1:fftsize/2-1]/fftsize*fs, imag(gamma_noise_hat), 'm--')
hold on; plot([1:fftsize/2-1]/fftsize*fs, gamma_noise, 'k-', 'LineWidth', 2)
legend('Re(est. coherence)', 'Im(est. coherence)', 'Theoretical coherence')
xlabel('Frequency (Hz)'); ylabel('Coherence')

number_frames = [4, 8, 16, 32, 64]
num_trials = 20000
Methods = {'GCC-PHAT', 'GCC-ML', 'MUSIC', 'Enhanced GCC-PHAT', 'GCC-PHAT, estimated noise spatial correlations', 'GEVD-MUSIC, estimated noise PSD matrices'}
Failed = zeros(length(Methods), length(number_frames));
for trial = 1 : num_trials
    for i_num_frames = 1 : length(number_frames)
        
        switch ceil(8*rand) % eight sentences, check the wav file; get some random continuous pieces for testing
            case 1
                t = round(202133 + rand*40453);
            case 2
                t = round(266858 + rand*33584);
            case 3
                t = round(326088 + rand*32134);
            case 4
                t = round(382952 + rand*39384);
            case 5
                t = round(444776 + rand*39843);
            case 6
                t = round(508128 + rand*28470);
            case 7
                t = round(562702 + rand*31981);
            case 8
                t = round(619260 + rand*38927);
        end
        
        rx1 = zeros(1, fftsize/2-1);
        rx2 = zeros(1, fftsize/2-1);
        cx12 = zeros(1, fftsize/2-1);
        for fn = 1 : number_frames(i_num_frames)
            X = fft(win.*x(:, t:t+fftsize-1), [], 2);
            rx1 = rx1 + X(1, 2:fftsize/2).*conj(X(1, 2:fftsize/2)); % the first bin is DC, thus discarded
            rx2 = rx2 + X(2, 2:fftsize/2).*conj(X(2, 2:fftsize/2));
            cx12 = cx12 + X(1, 2:fftsize/2).*conj(X(2, 2:fftsize/2));
            t = t + 160;
        end
        rx1=rx1/number_frames(i_num_frames); rx2=rx2/number_frames(i_num_frames); cx12=cx12/number_frames(i_num_frames);
         
        
        %% GCC-PHAT
        r = cx12./abs(cx12);
        r = ifftshift(ifft([0, r, 0, conj(r(end:-1:1))]));
        r = r(fftsize/2+1-delay_search_range : fftsize/2+1+delay_search_range);
        [~, tau_hat] = max(r);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(1, i_num_frames) = Failed(1, i_num_frames) + 1;
        end
        
        %% GCC-ML
        abs_gamma_hat2 = (abs(cx12)./sqrt(rx1.*rx2)).^2;
        r = cx12./(abs(cx12)).*abs_gamma_hat2./(1 - abs_gamma_hat2);
        r = ifftshift(ifft([0, r, 0, conj(r(end:-1:1))]));
        r = r(fftsize/2+1-delay_search_range : fftsize/2+1+delay_search_range);
        [~, tau_hat] = max(r);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(2, i_num_frames) = Failed(2, i_num_frames) + 1;
        end
        
        %% MUSIC. It can be slow due to solving lots of EVDs. You may comment out this part to skip it
        Music = zeros(1, 2*delay_search_range + 1);
        for k = 1 : fftsize/2-1
            P = [rx1(k), cx12(k); conj(cx12(k)), rx2(k)];
            [V,D] = eig(P);
            d = real(diag(D)); % eigenvalues should be nonnegative, although P is not Hermitian here
            if d(1) > d(2)
                v = V(:,2);
            else
                v = V(:,1);
            end
            for tau = -delay_search_range : delay_search_range
                a = [1; exp(sqrt(-1)*2*pi*k/fftsize*tau)];
                Music(tau+delay_search_range+1) = Music(tau+delay_search_range+1) + (abs(a'*v))^2; % some use abs(a'*v). I follow the standard MUSIC implementation
            end
        end
        Music = (fftsize/2-1)./Music;
        [~, tau_hat] = max(Music);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(3, i_num_frames) = Failed(3, i_num_frames) + 1;
        end
        
        %% enhanced PHAT
        gamma_hat = cx12./sqrt(rx1.*rx2);
        lambda = (1.0 - real(gamma_noise.*conj(gamma_hat)))./(1.0 - gamma_noise.*conj(gamma_noise));
        lambda = lambda - sqrt(lambda.^2 - (1.0 - gamma_hat.*conj(gamma_hat))./(1.0 - gamma_noise.*conj(gamma_noise))); % Pv in the paper
        gamma_ml = (gamma_hat - lambda.*gamma_noise)./(1 - lambda);
        r = ifftshift(ifft([0, gamma_ml, 0, conj(gamma_ml(end:-1:1))]));
        r = r(fftsize/2+1-delay_search_range : fftsize/2+1+delay_search_range);
        [~, tau_hat] = max(r);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(4, i_num_frames) = Failed(4, i_num_frames) + 1;
        end
        
        %% GCC-PHAT with estimated noise statistics
        r = cx12 - cn12;
        r = r./abs(r);
        r = ifftshift(ifft([0, r, 0, conj(r(end:-1:1))]));
        r = r(fftsize/2+1-delay_search_range : fftsize/2+1+delay_search_range);
        [~, tau_hat] = max(r);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(5, i_num_frames) = Failed(5, i_num_frames) + 1;
        end
        
        %% GEVD-MUSIC. It can be slow due to solving lots of EVDs. You may comment out this part to skip it
        Music = zeros(1, 2*delay_search_range + 1);
        for k = 1 : fftsize/2-1
            P = inv( [rn1(k), cn12(k); conj(cn12(k)), rn2(k)] ) ...
                   * [rx1(k), cx12(k); conj(cx12(k)), rx2(k)];
            [V,D] = eig(P);
            d = real(diag(D)); % eigenvalues should be nonnegative, although P is not Hermitian here
            if d(1) > d(2)
                v = V(:,2);
            else
                v = V(:,1);
            end
            for tau = -delay_search_range : delay_search_range
                a = [1; exp(sqrt(-1)*2*pi*k/fftsize*tau)];
                Music(tau+delay_search_range+1) = Music(tau+delay_search_range+1) + (abs(a'*v))^2; % some use abs(a'*v). I follow the standard MUSIC implementation
            end
        end
        Music = (fftsize/2-1)./Music;
        [~, tau_hat] = max(Music);
        if abs(tau_hat - delay_search_range-1-tau_true) >= 1
            Failed(6, i_num_frames) = Failed(6, i_num_frames) + 1;
        end
    end
    Failed
end
figure;
semilogy(10*number_frames, 100*Failed(1,:)/num_trials, 'p-b')
hold on; semilogy(10*number_frames, 100*Failed(2,:)/num_trials, 'x-m')
hold on; semilogy(10*number_frames, 100*Failed(3,:)/num_trials, '*-r')
hold on; semilogy(10*number_frames, 100*Failed(4,:)/num_trials, '+-k')
hold on; semilogy(10*number_frames, 100*Failed(5,:)/num_trials, 's-c')
hold on; semilogy(10*number_frames, 100*Failed(6,:)/num_trials, 'o-g')
xlabel('Length of speech signal (ms)')
ylabel('Percentage of failed TDOA estimates')
legend('GCC-PHAT', 'GCC-ML', 'MUSIC', 'Enhanced GCC-PHAT', 'GCC-PHAT, estimated noise spatial correlations', 'GEVD-MUSIC, estimated noise PSD matrices')