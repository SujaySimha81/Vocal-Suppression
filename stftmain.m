[y,fs]=audioread('Recording 5.wav');
y2=y(1:14539,1);
y1=transpose(y2);
t=size(y1);
s=t(1,2);
x=y1(1,1:s);

M=801;
w1=hamming(M);
w=transpose(w1);
N=1024;
H=40;

[mx,px]=stftanalysis(x,w,N,H);

