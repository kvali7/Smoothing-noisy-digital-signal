%snr  cables

clear all
% clearvars -except ipts1 ipts3 ipts4 iptl1 iptl3 iptl4;
clc;

%channels 9, 10, 11, 12
% columna 9, 10, 11, 12
%SPO2
%PLETH
%ETCO2
%PULSE
%no header
% no time
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
A = importdata('Test Run Data 1.txt');

startA = 100000;

N = size(A(startA:end,:),1);

f = 0:Fs/N:Fs-1/N;
t = 0:deltat:deltat*N-deltat;


figure (1)

for i=1:4
    
   subplot(4,1,i)
   
   plot(t, A(startA:end,8+i))
   
   title(cell2mat(titles1(i)));
   xlabel('Time (s)')
   ylabel ('Volts (V)')
end

figure(3);
for i=2
    subplot(4,1,i)
    plot(f, mag2db(abs(fft(hanning(N).*A(startA:end,8+i)/N))));
    xlim([0 Fs/2])
    xlabel('Frequency (Hz)')
end


[X1, X3, X4] = estimate_discrete_signal (threshold, A(startA:end,9), A(startA:end,11), A(startA:end,12));


% if ~(exist('ipts1')==1)
%     [ipts1,~] = findchangepts(A(startA:end,9),'MaxNumChanges',6,'Statistic','mean');
%     ipts1 = ipts1 -1;
%     ipts1 = [ipts1; size(A(startA:end,9),1)];
% end
% if ~(exist('ipts3')==1)
%     [ipts3,~] = findchangepts(A(startA:end,11),'MaxNumChanges',12,'Statistic','mean');
%     ipts3 = ipts3 -1;
%     ipts3 = [ipts3; size(A(startA:end,11),1)];
% end
% if ~(exist('ipts4')==1)
%     [ipts4,~] = findchangepts(A(startA:end,12),'MaxNumChanges',50,'Statistic','mean');
%     ipts4 = ipts4 -1;
%     ipts4 = [ipts4; size(A(startA:end,12),1)];
% end


% X1=[];
% temp=startA;
% for j= 1:7
%     X1 = [X1; mean(A(temp:startA+ipts1(j)-1,9))*ones(startA+ipts1(j)-temp,1)];
%     temp = startA + ipts1(j);
% 
% end
% 
% % X1=[];
% % temp=start;
% % X1 = [X1; mean(A(temp:start+4727-1,9))*ones(start+4727-temp,1)];
% % temp = start + 4727;
% % 
% % X1 = [X1; mean(A(temp:start+6227-1,9))*ones(start+6227-temp,1)];
% % temp = start + 6227;
% % 
% % X1 = [X1; mean(A(temp:start+10227-1,9))*ones(start+10227-temp,1)];
% % temp = start + 10227;
% % 
% % X1 = [X1; mean(A(temp:start+16227-1,9))*ones(start+16227-temp,1)];
% % temp = start + 16227;
% % 
% % X1 = [X1; mean(A(temp:start+27226-1,9))*ones(start+27226-temp,1)];
% % temp = start + 27226;
% % 
% % X1 = [X1; mean(A(temp:start+30726-1,9))*ones(start+30726-temp,1)];
% % temp = start + 30726;
% % 
% % X1 = [X1; mean(A(temp:start+50350-1,9))*ones(start+50350-temp,1)];
% % temp = start + 50350;
% 
% 
% X3=[];
% temp=startA;
% for j= 1:13
%     X3 = [X3; mean(A(temp:startA+ipts3(j)-1,11))*ones(startA+ipts3(j)-temp,1)];
%     temp = startA + ipts3(j);
% 
% end
% 
% 
% X4=[];
% temp=startA;
% for j= 1:51
%     X4 = [X4; mean(A(temp:startA+ipts4(j)-1,12))*ones(startA+ipts4(j)-temp,1)];
%     temp = startA + ipts4(j);
% 
% end



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

maeshortch1 = mae(X1, A(startA:end,9))
maeshortch3 = mae(X3, A(startA:end,11))
maeshortch4 = mae(X4, A(startA:end,12))

mseshortch1 = mse(X1, A(startA:end,9))
mseshortch3 = mse(X3, A(startA:end,11))
mseshortch4 = mse(X4, A(startA:end,12))

snrshortch1 = snr(X1, A(startA:end,9)-X1)
snrshortch3 = snr(X3, A(startA:end,11)-X3)
snrshortch4 = snr(X4, A(startA:end,12)-X4)

maxeshortch1 = max(abs(A(startA:end,9)-X1))
maxeshortch3 = max(abs(A(startA:end,11)-X3))
maxeshortch4 = max(abs(A(startA:end,12)-X4))


%% long cable
B = importdata('Test Run Data 2.txt');
startB = 1;

N = size(B(startB:end,:),1);

f = 0:Fs/N:Fs-1/N;
t = 0:deltat:deltat*N-deltat;

figure (2)

for i=1:4
    
   subplot(4,1,i)
   
   plot(t, B(startB:end,8+i))
   
    title(cell2mat(titles2(i)));
    xlabel('Time (s)')
    ylabel ('Volts (V)')
end

figure(4);
for i=2
    subplot(4,1,i)
    plot(f, mag2db(abs(hanning(N).*fft(B(startB:end,8+i)/N))));
    xlim([0 Fs/2])
    xlabel('Frequency (Hz)')
end


[Y1, Y3, Y4] = estimate_discrete_signal (threshold, B(startB:end,9), B(startB:end,11), B(startB:end,12));


% if ~(exist('iptl1')==1)
%     [iptl1,~] = findchangepts(B(startB:end,9),'MaxNumChanges',4,'Statistic','mean');
%     iptl1 = iptl1 -1;
%     iptl1 = [iptl1; size(B(startB:end,9),1)];
% end
% if ~(exist('iptl3')==1)
%     [iptl3,~] = findchangepts(B(startB:end,11),'MaxNumChanges',20,'Statistic','mean');
%     iptl3 = iptl3 -1;
%     iptl3 = [iptl3; size(B(startB:end,11),1)];
% end
% if ~(exist('iptl4')==1)
%     [iptl4,~] = findchangepts(B(startB:end,12),'MaxNumChanges',54,'Statistic','mean');
%     iptl4 = iptl4 -1;
%     iptl4 = [iptl4; size(B(startB:end,12),1)];
% end
% 
% 
% Y1=[];
% temp=startB;
% for j=1:5
%     Y1 = [Y1; mean(B(temp:startB+iptl1(j)-1,9))*ones(startB+iptl1(j)-temp,1)];
%     temp = startB + iptl1(j);
% end
% 
% 
% % Y1=[];
% % temp=startB;
% % Y1 = [Y1; mean(B(temp:startB+5371-1,9))*ones(startB+5371-temp,1)];
% % temp = startB + 5371;
% % 
% % Y1 = [Y1; mean(B(temp:startB+7588-1,9))*ones(startB+7588-temp,1)];
% % temp = startB + 7588;
% % 
% % Y1 = [Y1; mean(B(temp:startB+48586-1,9))*ones(startB+48586-temp,1)];
% % temp = startB + 48586;
% % 
% % Y1 = [Y1; mean(B(temp:startB+52585-1,9))*ones(startB+52585-temp,1)];
% % temp = startB + 52585;
% % 
% % Y1 = [Y1; mean(B(temp:startB+65828-1,9))*ones(startB+65828-temp,1)];
% % temp = startB + 65828;
% 
% 
% Y3=[];
% temp=startB;
% for j=1:21
%     Y3 = [Y3; mean(B(temp:startB+iptl3(j)-1,11))*ones(startB+iptl3(j)-temp,1)];
%     temp = startB + iptl3(j);
% end
% 
% 
% Y4=[];
% temp=startB;
% for j=1:55
%     Y4 = [Y4; mean(B(temp:startB+iptl4(j)-1,12))*ones(startB+iptl4(j)-temp,1)];
%     temp = startB + iptl4(j);
% end


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

maelongch1 = mae(Y1, B(startB:end,9))
maelongch3 = mae(Y3, B(startB:end,11))
maelongch4 = mae(Y4, B(startB:end,12))

mselongch1 = mse(Y1, B(startB:end,9))
mselongch3 = mse(Y3, B(startB:end,11))
mselongch4 = mse(Y4, B(startB:end,12))

snrlongch1 = snr(Y1, B(startB:end,9)-Y1)
snrlongch3 = snr(Y3, B(startB:end,11)-Y3)
snrlongch4 = snr(Y4, B(startB:end,12)-Y4)

maxelongch1 = max(abs(B(startB:end,9)-Y1))
maxelongch3 = max(abs(B(startB:end,11)-Y3))
maxelongch4 = max(abs(B(startB:end,12)-Y4))


% AACQ = load_acq('Test Run Data 1.acq')
