function angle = srp_phat(filename, resolution)
%% This code is made by Seongbin kim 19-05-01
%% for estimate angle of 2ch recorded audio signal
%% angle is -90 to 90 degrees that left 90 degree is -90, right 90 is +90
%% Inha university Bio IT Lab
%% email : zelabean@naver.com
[signal, fs] = audioread(filename);
%% initial params
w_size = 25/1000;                                       %% window size : 25 ms
shift_size = 10/1000;                                   %% frame shift size : 10 ms
stft_point = 256;                                       %% stft point : 256
d = 0.094;                                              %% distance of 2 microphone (m)
c = 340;                                                %% wind speed (m/s)
pre_angle = -90:resolution:90;                          %% pre_angle is created angle by resolution, ex) resolution 1 -> -90, -89, -88...; resolution 0.1 -> -90, -89.9 ...
tau = d*sin(deg2rad(pre_angle))/c;
rightSig = signal(:,1);                                 %% split 2ch signal to 1ch signal
leftSig = signal(:,2);
X1 = spectrogram(rightSig,hamming(w_size*fs),shift_size*fs,stft_point,fs, 'yaxis');
X2 = spectrogram(leftSig,hamming(w_size*fs),shift_size*fs,stft_point,fs, 'yaxis');
G = X1.*conj(X2);
psi = 1./abs(G);

tmp = size(X1);
n = tmp(2);

p = length(tau);
R = zeros(p,n);
for i = 1:p
    sum = 0;
    for k = 1:129
        sum = sum + psi(k,:).*G(k,:)*exp(1j*2*pi*(fs/stft_point)*k*tau(i));
    end
    R(i,:) = sum;
end
[a,b] = max(R); %% find argmax
angle = mean(rad2deg(asin(tau(b)/d*c))); %% estimate angle

end
