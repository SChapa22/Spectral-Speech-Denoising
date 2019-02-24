% Audio Denoising by Time-Frequency Block Thresholding
% Guoshen Yu, Stéphane Mallat, Fellow, IEEE, and Emmanuel Bacry%
% IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 56, NO. 5, MAY 2008

%-----Completed on 20.12.2018 - 20:20
%-----Coded by: Y.Emir ARITÜRK
%-----Version No: 1218-2-2

%   load noisy sound example
fprintf('-> Step 1/5: Load noisy.wav:\n');
[y,Fs]=audioread('your path to input sound source'); % y=sampled data, Fs=sampling frequency
x=y(100300:end);  %x is the "beginning removed" and transposed version of the sampled data y. 
fprintf('Load noisy.wav: OK...\n');


%   STFT parameters for spectogram (check explanation of spectogram function
NFFT=1024;
window_length=round(0.03*Fs); %window length for stft
window=kaiser(window_length, 15); %create a kaiser window in length of window_length. the input signal x will be divided into segments in length of window variable.
overlap=floor(0.45*window_length); %number of window samples without overlapping

%   Signal parameters
t_min=0.001;    %the moment getting the noise for the first time (sn)
t_max=1.1;   % noise at the last moment(sn)

%   construct spectrogram
[S,F,T] = spectrogram(x, window, window_length-overlap, NFFT,Fs); %S: short time fourier transform, NFFT: number of DFT points to be used in STFT
[Nf,Nw]=size(S); %the frequency and time components of STFT; F:frequency vector and T:time vector has been scalarized into Nf and Nw respectively


%   noisy spectrum extraction
fprintf('-> Step 2/5: Extract noise spectrum \n');
t_index= T>t_min & T<t_max; %inspection interval
absS_noise=abs(S(:,t_index)).^2; %power spectrum of the noise
noise_spectrum=mean(absS_noise,2); %the resulting “denoised” signal estimator. it'll be used later in SNR estimation and linked to Wiener attenuation
noise_specgram=repmat(noise_spectrum,1,Nw); %creating a (repeat copy) noise matrix from noise spectrum vector to make size same as S. This value will be used in SNR estimation.
fprintf('Extract noise spectrum: OK...\n');


%   parameters for diagonal SNR Estimation
fprintf('-> Step 3/5: Estimate SNR -');
SNR_estimator=max(((abs(S).^2)./noise_specgram)-1,0); %(abs(S).^2)./noise_specgram)is a posteriori SNR defined by a priori SNR (abs(S).^2). 
% SNR_estimator will help to find the Wiener attenuation in following parts
fprintf('Estimate SNR: OK...\n');


%   Parameters for Wiener attenuation rule 
beta1=0.5; %Must be greater than or equal to 0
beta2=1;  %Must be greater than or equal to 0
lambda=3; %Must be greater than 1. Lambda is an over-subtraction factor to compensate variation of noise amplitude. 

%   Compute time-frequency attenuation map
fprintf('-> Step 4/5: Compute TF attenuation map -\n');
attenuation_factor=max((1-lambda*((1./(SNR_estimator+1)).^beta1)).^beta2,0); %Wiener attenuation rule
temp_S=attenuation_factor.*S; %temp S value will be used in inverse STFT when getting the speech-only signal
fprintf('Compute TF attenuation map: OK...\n');


%   Compute Inverse STFT   
fprintf('-> Step 5/5: Compute Inverse STFT:\n');
ind=mod((1:window_length)-1,Nf)+1;
output_signal=zeros((Nw-1)*overlap+window_length,1);

for indice=1:Nw %Overlap add technique to get the speech only signal
    left_index=((indice-1)*overlap) ;
    index=left_index+[1:window_length];
    temp_ifft=real(ifft(temp_S(:,indice),NFFT)); %inverse fast-fourier transform
    output_signal(index)= output_signal(index)+temp_ifft(ind).*window;
end
fprintf('Compute Inverse STFT: OK...\n');

%--------------------------------DISPLAY----------------------------------

%Display the the original signal
figure
subplot(2,1,1);
plot ((1:length(x))/Fs,x)
xlabel('Time (s)');
ylabel('Amplitude');
t_index=find(T>t_min & T<t_max);
plot((1:length(x))/Fs,x);
xlabel('Time (s)');
ylabel('Amplitude');
hold on;

%show noise in red color
noise_interval=(T(t_index(1))*Fs : T(t_index(end))*Fs);
plot(noise_interval/Fs,x(noise_interval),'r');
legend('Original signal','Audience noise');
title('Original Sound');
hold off;

%    show denoised signal
subplot(2,1,2);
plot((1:length(output_signal))/Fs,output_signal);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Denoised signal');
title('Sound without audience noise');

%show spectrogram
t_epsilon=0.001;
figure
S_one_sided=max(S(1:length(F)/2,:),t_epsilon); %keep only the positive frequency
pcolor(T,F(1:end/2),10*log10(abs(S_one_sided))); 
shading interp;
colormap('jet');
title('Kaiser Windowed Spectrogram Spectrogram (\beta = 15): speech + audience');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure
S_one_sided=max(temp_S(1:length(F)/2,:),t_epsilon); %keep only the positive frequency
pcolor(T,F(1:end/2),10*log10(abs(S_one_sided))); 
shading interp;
colormap('jet');
title('Kaiser Windowed Spectrogram (\beta = 15)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%---- Result Speeches ---
wait_time = 5; %total time(s) for waiting from start of the first sound file

%Play original sound
fprintf('\nPlaying noisy speech:\n');
original_speech = audioplayer(x(1:(wait_time-1)*Fs), Fs);
play(original_speech);
fprintf('Noisy speech play completed...\n');

pause(wait_time)

%Play denoised sound
fprintf('\nPlaying denoised speech:');
denoised_speech = audioplayer(output_signal(1:(wait_time-1)*Fs), Fs);
play(denoised_speech,Fs);
fprintf('Denoised speech play completed...\n');
audiowrite('your path to output sound source',output_signal,Fs);
fprintf('denoised.wav file has been successfuly created...\n');
