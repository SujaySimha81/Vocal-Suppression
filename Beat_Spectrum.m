clc
close all
tic

% load a .wav file
[x, fs] = audioread('Mono_Humma.wav'); % get the samples of the .wav file
a = size(x);
%x = (x(:, 1) + x(:, 2))/2;                        % get the first channel
x = x(:,1);
xmax = max(abs(x));                 % find the maximum abs value
x = x/xmax;                         % scaling the signal

% define analysis parameters
xlen = length(x);                   % length of the signal
wlen = 1024;                        % window length (recomended to be power of 2)
h = wlen/2;                         % hop size (recomended to be power of 2)
nfft = 2048;                        % number of fft points (recomended to be power of 2)

% define the coherent amplification of the window
% K = sum(hamming(wlen, 'periodic'))/wlen;

% perform STFT
[s, f, t] = STFT_1(x, wlen, h, nfft, fs);
dim = size(s);
n_rows = dim(1);
n_cols = dim(2);

V = (abs(s)).^2;
B = zeros(n_rows , n_cols);
toc
tic
% for i = 1:n_rows
%     for j = 1:n_cols
%         for k = 1:n_cols - j + 1
%             B(i,j) = B(i,j) + V(i,k).*V(i,k+j-1);
%         end
%         B(i,j) = B(i,j).*(1/(n_cols - j + 1));
%     end
% end
for j = 1:n_cols
    for k = 1:n_cols - j + 1
        B(:,j) = B(:,j) + V(:,k).*V(:,k+j-1);
    end
    %B(:,j) = sum((V(:,1:n_cols-j+1).*V(:,j:n_cols)),2);
    B(:,j) = B(:,j).*(1/(n_cols - j + 1));
end
toc
tic
b = zeros(1,n_cols);
end_of_song = max(t);
del = ceil((0.5*n_cols)/end_of_song);
b(:) = sum(B);
b = b(del:length(b));
b = b./max(b);
temp = floor(length(b)/4);
l = length(b) - temp;
b = b(1:l);
[maximum , time_max] = max(b);
%b = b - min(b);
time = t(del:length(t));
time = time(1:l);

figure
plot(time,b);

%% Finding repeating period p
test = b(1:time_max - 1);
[t,loc] = max(test);
if (maximum - t) > 0.1
    min_height = maximum - 0.16;
else
    min_height = maximum - 0.1;
end 
[local_peaks,locations] = findpeaks(b(1:ceil((time_max)/2)),'MinPeakHeight',min_height);
peak_time = time(time_max);
time_locations = time(locations);
temp = zeros(1,length(locations));
temp1 = zeros(1,length(locations));
check = zeros(1,3);
for i = 1:length(time_locations)
%     plus_minus_correction = [time_locations(i) - 0.015, time_locations(i), time_locations(i) + 0.015];
%     for j = 1:3
%         check(j) = abs(peak_time/plus_minus_correction(j) - round(peak_time/plus_minus_correction(j)));
%     end
%     [value,index] = min(check);
%    temp(i) = check(index);
    temp(i) = abs(peak_time/time_locations(i) - round(peak_time/time_locations(i)));
%    temp1(i) = plus_minus_correction(index);
    temp(i) = (0.97*temp(i) + 0.03*time_locations(i));
end
[dummy , p] = min(temp);
if (isempty(temp) == false) % Change made here
    period = time_locations(p);
else
    period = peak_time;
end
%period = time_locations(p);
toc
