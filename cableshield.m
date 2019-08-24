%snr  cables

clear all
% clearvars -except ipts1 ipts3 ipts4 iptl1 iptl3 iptl4;
clc;

%channels 13, 14, 15, 16
%columns 11, 12, 13, 14
%SPO2
%PLETH
%ETCO2
%PULSE

%with header
% column 1 is time in minutes
Fs = 500;
deltat = 1/ Fs;

threshold = 0.0001;




titles1 = {'SPO2 short cable';...
    'PLETH short cable';...
    'ETCO2 short cable';...
    'PULSE short cable'};

titles2 = {'SPO2 long cable';...
    'PLETH long cable';...
    'ETCO2 long cable';...
    'PULSE long cable'};


titlesest1 = {'SPO2 short cable estimate';...
    'PLETH short cable estimate';...
    'ETCO2 short cable estimate';...
    'PULSE short cable estimate'};

titlesest2 = {'SPO2 long cable estimate';...
    'PLETH long cable estimate';...
    'ETCO2 long cable estimate';...
    'PULSE long cable estimate'};

%% short cable
A = importdata('8-22-19 test run 1 short cable.txt');
Aheaders = A.colheaders;
disp (Aheaders);
A = A.data (2:end,:);
startA = 36397;

N = size(A(startA:end,:),1);

f = 0:Fs/N:Fs-1/N;
t = 0:deltat:deltat*N-deltat;

figure (1)

for i=1:4
    
   subplot(4,1,i)
   
   plot(t, A(startA:end,10+i))
   
   title(cell2mat(titles1(i)));
   xlabel('Time (s)')
   ylabel ('Volts (V)')
end

figure(3);
for i=2
    subplot(4,1,i)
    plot(f, mag2db(abs(fft(hanning(N).*A(startA:end,10+i)/N))));
    xlim([0 Fs/2])
    xlabel('Frequency (Hz)')
end

numchangepoints1 = 13;
numchangepoints3 = 8;
numchangepoints4 = 55;

if ~(exist('ipts1')==1)
    [ipts1,~] = findchangepts(A(startA:end,11),'MinThreshold',threshold,'Statistic','mean');
    ipts1 = ipts1 -1;
    ipts1 = [ipts1; size(A(startA:end,11),1)];
end
if ~(exist('ipts3')==1)
    [ipts3,~] = findchangepts(A(startA:end,13),'MinThreshold',threshold,'Statistic','mean');
    ipts3 = ipts3 -1;
    ipts3 = [ipts3; size(A(startA:end,13),1)];
end
if ~(exist('ipts4')==1)
    [ipts4,~] = findchangepts(A(startA:end,14),'MinThreshold',threshold,'Statistic','mean');
    ipts4 = ipts4 -1;
    ipts4 = [ipts4; size(A(startA:end,14),1)];
end




X1=[];
temp=startA;
for j= 1:numchangepoints1 + 1
    X1 = [X1; mean(A(temp:startA+ipts1(j)-1,11))*ones(startA+ipts1(j)-temp,1)];
    temp = startA + ipts1(j);

end


X3=[];
temp=startA;
for j= 1:numchangepoints3 + 1
    X3 = [X3; mean(A(temp:startA+ipts3(j)-1,13))*ones(startA+ipts3(j)-temp,1)];
    temp = startA + ipts3(j);

end


X4=[];
temp=startA;
for j= 1:numchangepoints4 + 1
    X4 = [X4; mean(A(temp:startA+ipts4(j)-1,14))*ones(startA+ipts4(j)-temp,1)];
    temp = startA + ipts4(j);

end

figure(5)

subplot(4,1,1)
plot(t, X1);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest1(1)));
subplot(4,1,3)
plot(t, X3);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest1(3)));
subplot(4,1,4)
plot(t, X4);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest1(4)));

maeshortch1 = mae(X1, A(startA:end,11))
maeshortch3 = mae(X3, A(startA:end,13))
maeshortch4 = mae(X4, A(startA:end,14))

mseshortch1 = mse(X1, A(startA:end,11))
mseshortch3 = mse(X3, A(startA:end,13))
mseshortch4 = mse(X4, A(startA:end,14))

snrshortch1 = snr(X1, A(startA:end,11)-X1)
snrshortch3 = snr(X3, A(startA:end,13)-X3)
snrshortch4 = snr(X4, A(startA:end,14)-X4)

maxeshortch1 = max(abs(A(startA:end,11)-X1))
maxeshortch3 = max(abs(A(startA:end,13)-X3))
maxeshortch4 = max(abs(A(startA:end,14)-X4))


[ests1, ests3, ests4] = estimate_discrete_signal (threshold, A(startA:end,11), A(startA:end,13), A(startA:end,14));

%% long cable
B = importdata('8-22-19 test run 1 long cable.txt');
Bheaders = B.colheaders;
disp (Bheaders);
B = B.data (2:end,:);

startB = 1;

N = size(B(startB:end,:),1);

f = 0:Fs/N:Fs-1/N;
t = 0:deltat:deltat*N-deltat;

figure (2)

for i=1:4
    
   subplot(4,1,i)
   
   plot(t, B(startB:end,10+i))
   
    title(cell2mat(titles2(i)));
    xlabel('Time (s)')
    ylabel ('Volts (V)')
end

figure(4);
for i=2
    subplot(4,1,i)
    plot(f, mag2db(abs(hanning(N).*fft(B(startB:end,10+i)/N))));
    xlim([0 Fs/2])
    xlabel('Frequency (Hz)')
end

numchangepointl1 = 9;
numchangepointl3 = 21;
numchangepointl4 = 68;

if ~(exist('iptl1')==1)
    [iptl1,~] = findchangepts(B(startB:end,11),'MinThreshold',threshold,'Statistic','mean');
    iptl1 = iptl1 -1;
    iptl1 = [iptl1; size(B(startB:end,11),1)];
end
if ~(exist('iptl3')==1)
    [iptl3,~] = findchangepts(B(startB:end,13),'MinThreshold',threshold,'Statistic','mean');
    iptl3 = iptl3 -1;
    iptl3 = [iptl3; size(B(startB:end,13),1)];
end
if ~(exist('iptl4')==1)
    [iptl4,~] = findchangepts(B(startB:end,14),'MinThreshold',threshold,'Statistic','mean');
    iptl4 = iptl4 -1;
    iptl4 = [iptl4; size(B(startB:end,14),1)];
end


Y1=[];
temp=startB;
for j=1:numchangepointl1 + 1
    Y1 = [Y1; mean(B(temp:startB+iptl1(j)-1,11))*ones(startB+iptl1(j)-temp,1)];
    temp = startB + iptl1(j);
end



Y3=[];
temp=startB;
for j=1:numchangepointl3 + 1 
    Y3 = [Y3; mean(B(temp:startB+iptl3(j)-1,13))*ones(startB+iptl3(j)-temp,1)];
    temp = startB + iptl3(j);
end


Y4=[];
temp=startB;
for j=1:numchangepointl4 + 1
    Y4 = [Y4; mean(B(temp:startB+iptl4(j)-1,14))*ones(startB+iptl4(j)-temp,1)];
    temp = startB + iptl4(j);
end


figure(6)

subplot(4,1,1)
plot(t, Y1);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest2(1)));
subplot(4,1,3)
plot(t, Y3);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest2(3)));
subplot(4,1,4)
plot(t, Y4);
   xlabel('Time (s)')
   ylabel ('Volts (V)')
title(cell2mat(titlesest2(4)));

maelongch1 = mae(Y1, B(startB:end,11))
maelongch3 = mae(Y3, B(startB:end,13))
maelongch4 = mae(Y4, B(startB:end,14))

mselongch1 = mse(Y1, B(startB:end,11))
mselongch3 = mse(Y3, B(startB:end,13))
mselongch4 = mse(Y4, B(startB:end,14))

snrlongch1 = snr(Y1, B(startB:end,11)-Y1)
snrlongch3 = snr(Y3, B(startB:end,13)-Y3)
snrlongch4 = snr(Y4, B(startB:end,14)-Y4)

maxelongch1 = max(abs(B(startB:end,11)-Y1))
maxelongch3 = max(abs(B(startB:end,13)-Y3))
maxelongch4 = max(abs(B(startB:end,14)-Y4))

[estl1, estl3, estl4] = estimate_discrete_signal (threshold, B(startB:end,11), B(startB:end,13), B(startB:end,14));
% AACQ = load_acq('Test Run Data 1.acq')
