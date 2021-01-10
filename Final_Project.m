% CSE 3313 - Final Project
% Nikolas Murguia, 1001666001, December 7,2020
clear all;
% Part 1: Audio Frequency Analysis

% a: Get sampling frequency of audio file
[y, Fs] = audioread('noisyaudio.wav');

% b: Find DFT of audio sample
X = fft(y);

% c: Form frequency axis w span -Fs/2 to Fs/2
Axis = linspace(-Fs/2, Fs/2, length(y));

% d: Plot the magnitude of DTF vs frecuency
subplot(4,1,1);
plot(Axis, fftshift(abs(X)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude vs Frequency');

% e: Plot the normalized log plot of DFT
x = abs(X)/max(abs(X));
dB = 20.*log10(x);

subplot(4,1,2);
plot(Axis, fftshift(dB));
xlabel('Frequency (Hz)');
ylabel('dB');
title('Normalized Log Plot');


% Part 2: Filter Design

% f: Decide wp and ws measured in Hz
wp = 2100;
ws = 2500;

% g: Assume attenuation for passband (ATp) = -1dB
%    use Normalized Log Plot to find attenuation of stopband (ATs)
ATp = -1;
ATs = -50;

% h: Create Butterworth analog filter
%    Filter order = 0.5*log(ks/kp)/log(ws/wp)
%    ks = 1/(10^(ATs/20))^2 - 1
%    kp = 1/(10^(ATp/20))^2 - 1
ks = 1/((10.^(ATs/20)).^2)-1;
kp = 1/((10.^(ATp/20)).^2)-1;
N = (log10(ks/kp)/(2*log10(ws/wp)));
N = round(double(N))+1;

% i: Find Cutoff Frequency
COFreq = wp/(10^(kp)^(1/(2*N)));

% j: Calculate Ha and plot the logarithmic gain
Ha = 20*log10(sqrt(1./(1+(Axis/COFreq).^(2*N))));
subplot(4,1,3);
plot(Axis, Ha);
xlabel('Frequency (Hz)');
ylabel('dB');
title('Logarithmic Gain of Frequency Response');

% Part 3: Filter Implementation

% 3a: Use butter() command to design filter
%     Calculate Wn
Wn = COFreq/(Fs/2);
[b,a] = butter(N,Wn);

% 3b: Use filter() command and the a&b coefficients
filterAudio = filter(b,a,y);

% 3c: Calculate DFT and plot magnitude
X2 = fft(filterAudio);
subplot(4,1,4);
plot(Axis, fftshift(abs(X2)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude (Filtered) vs Frequency');

% 3e: Play the unfiltered and then filtered audio
sound(y,Fs);            % Unfiltered
pause(13);
sound(filterAudio, Fs); % Filtered

% 3f: Write the new filtered audio
audiowrite('filterdaudio.wav', filterAudio, Fs);



