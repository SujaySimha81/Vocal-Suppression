% load a .mp3 file
[x, fs] = audioread('song6.mp3'); % get the samples of the .mp3 file
a = size(x);
%x = (x(:, 1) + x(:, 2))/2;              
x = x(:,1);  % get the first channel
x = x(310144:2310144,1);
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
V=abs(s);
Vs = (abs(s)).^2;
B = zeros(n_rows , n_cols);

for j = 1:n_cols
    for k = 1:n_cols - j + 1
        B(:,j) = B(:,j) + Vs(:,k).*Vs(:,k+j-1);
    end
    
    B(:,j) = B(:,j).*(1/(n_cols - j + 1));
end

b = zeros(1,n_cols);
b(:)=sum(B);
b=b/n_rows;
% b=b/b(1,1);
b=b/max(b);
new=round(3*n_cols/4);
b1=zeros(1,new);
k=1;
for i=1:n_cols
    if(b(1,i)<0.75)
        b1(1,k)=b(1,i);
        k=k+1;
    end
end
end_of_song = max(t);

% del = ceil((0.5*n_cols)/end_of_song);
% b(:) = sum(B);
% b = b(del:length(b));
% b = b./max(b);
% temp = floor(length(b)/4);
% l = length(b) - temp;
% 
% b = b(1:l);

% [maximum , time_max] = max(b);
% %b = b - min(b);
% time = t(del:length(t));
% time = time(1:l);

figure
time=0:end_of_song/n_cols:end_of_song;
time=time(1:n_cols);
plot(time,b);

%% Finding repeating period p
l=new;
% del = ceil((0.5*n_cols)/end_of_song);
r=floor(l/3);
J=zeros(1,r);
for j=1:floor(n_cols/3)
    del=floor(j/4);
    delta=floor(3*j/4);
    I=0;
    for k=1:floor(0.2*n_cols/j)
        i=k*j;
        T=b1(1,i-del:i+del);
        T1=b1(1,i-delta:i+delta);
        [q1,h]=max(T);
        [q2,h1]=max(T1);
        if(h==h1)
            I=I+b1(h)-mean(T1);
        end
    end
    J(1,j)=I/floor(l/j);
end
[m_ele,p]=max(J);

%% repeating segment modelling
ps=round(p*n_cols/end_of_song);
%p=p;
r=round(n_cols/ps);
S=zeros(n_rows,ps);
for i=1:n_rows
    for l=1:ps
        %M1=zeros(1,r);
        k=1:r;
        M1=V(i,l+(k-1)*ps);
        med=median(M1);
        S(i,l)=med;
        
    end
end


%% extracting repeating patterns

W=zeros(n_rows,n_cols);
for i=1:n_rows
    for l=1:ps
        k=1:r;
        W(i,l+(k-1)*ps)=min(S(i,l),V(i,l+(k-1)*ps));
    end
end

%M=zeros(n_rows,n_cols);
M=W./V;
I=M.*s;

[res, t1] = ISTFT_1(I, wlen, h, nfft, fs);
        
        
















