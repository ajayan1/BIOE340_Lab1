%% Part 1
clear; clc; close all
% 1)
healthy = readtable("Lab2_Healthy_Data_ECG.xlsx");
disease = readtable("Lab2_Disease_Data_ECG.xlsx");

healthy_time = healthy.Time;
healthy_leadI = healthy.LeadI./1000;
healthy_leadII = healthy.LeadII./1000;
healthy_leadIII = healthy.LeadIII./1000;
healthy_aVR = healthy.aVR./1000;
healthy_aVL = healthy.aVL./1000;
healthy_aVF = healthy.aVF./1000;
subplot_titles_1 = ["Healthy Lead I","Healthy Lead II","Healthy Lead III","Healthy Lead aVR","Healthy Lead aVL","Healthy Lead aVF"];

healthy_leads = [healthy_leadI,healthy_leadII,healthy_leadIII,healthy_aVR,healthy_aVL,healthy_aVF];
detrend_healthy_leads = [];
for i = 1:6
    [p, ~, mu] = polyfit(healthy_time,healthy_leads(:,i),4);
    detrend_healthy_leads(:,end+1) = healthy_leads(:,i) - polyval(p, healthy_time, [], mu);
end

figure(Name = '6-Lead ECG Healthy')
for i = 1:6
    subplot(2,3,i)
    plot(healthy_time,detrend_healthy_leads(:,i))
    title(subplot_titles_1(i))
    xlabel('Time (s)')
    ylabel('mV')
end

% figure(1)
% subplot(6,1,1)
% plot(healthy_time,healthy_leadI)
% title('Healthy Lead I')
% 
% subplot(6,1,2)
% plot(healthy_time,healthy_leadII)
% title('Healthy Lead II')
% 
% subplot(6,1,3)
% plot(healthy_time,healthy_leadIII)
% title('Healthy Lead III')
% 
% subplot(6,1,4)
% plot(healthy_time,healthy_aVR)
% title('Healthy Lead aVR')
% 
% subplot(6,1,5)
% plot(healthy_time,healthy_aVL)
% title('Healthy Lead aVL')
% 
% subplot(6,1,6)
% plot(healthy_time,healthy_aVF)
% title('Healthy Lead aVF')

% 3)
Cal_LeadIII = detrend_healthy_leads(:,2) - detrend_healthy_leads(:,1);
Cal_aVF = ((2.*detrend_healthy_leads(:,2))-detrend_healthy_leads(:,1))./(sqrt(3));
Cal_aVL = ((2.*detrend_healthy_leads(:,1))-detrend_healthy_leads(:,2))./(sqrt(3));
Cal_aVR = -(detrend_healthy_leads(:,2)+detrend_healthy_leads(:,1))./(sqrt(3));

a = isequal(Cal_LeadIII,detrend_healthy_leads(:,3))

% 4)
figure(Name = 'Derived 6-Lead ECG')
compare_derived = [Cal_LeadIII,detrend_healthy_leads(:,3),Cal_aVR,detrend_healthy_leads(:,4),Cal_aVL,detrend_healthy_leads(:,5),Cal_aVF,detrend_healthy_leads(:,6)];
subplot_titles_2 = ["Derived Lead III","Healthy Lead III","Derived Lead aVR","Healthy Lead aVR","Derived Lead aVL","Healthy Lead aVL","Derived Lead aVF","Healthy Lead aVF"];
for i = 1:8
    subplot(4,2,i)
    plot(healthy_time,compare_derived(:,i))
    title(subplot_titles_2(i))
    xlabel('Time (s)')
    ylabel('mV')
end
for i = 1:2:8
    b = isequal(compare_derived(:,i),compare_derived(:,i+1))
end

% subplot(4,2,1)
% plot(healthy_time,Cal_LeadIII)
% title('Derived Lead III')
% 
% subplot(4,2,2)
% plot(healthy_time,detrend_healthy_leads(:,3))
% title('Healthy Lead III')
% 
% subplot(4,2,3)
% plot(healthy_time,Cal_aVR)
% title('Derived Lead aVR')
% 
% subplot(4,2,4)
% plot(healthy_time,healthy_aVR)
% title('Healthy Lead aVR')
% 
% subplot(4,2,5)
% plot(healthy_time,Cal_aVF)
% title('Derived Lead aVF')
% 
% subplot(4,2,6)
% plot(healthy_time,healthy_aVF)
% title('Healthy Lead aVF')
% 
% subplot(4,2,7)
% plot(healthy_time,Cal_aVL)
% title('Derived Lead aVL')
% 
% subplot(4,2,8)
% plot(healthy_time,healthy_aVL)
% title('Healthy Lead aVL')

%% Part 2

% Average leads/detrend data
average_healthy_lead = mean(detrend_healthy_leads,2);
[p, ~, mu] = polyfit(healthy_time,average_healthy_lead,4);
detrend_avg_healthy = average_healthy_lead - polyval(p, healthy_time, [], mu);

% Initialize PQRST arrays
P_peaks = [];
P_locs = [];
Q_peaks = [];
Q_locs = [];
T_peaks = [];
T_locs = [];
S_peaks = [];
S_locs = [];
R_peaks = [];
R_locs = [];

% Find PRT peaks
[PRT_peaks,PRT_locs] = findpeaks(detrend_avg_healthy,NPeaks=25,MinPeakHeight=0.01,MinPeakDistance=20);

for i = 1:length(PRT_peaks)
    if mod(i-1,3) == 0
        T_peaks(end+1) = PRT_peaks(i);
        T_locs(end+1) = PRT_locs(i);
    elseif  mod(i-2,3) == 0
        P_peaks(end+1) = PRT_peaks(i);
        P_locs(end+1) = PRT_locs(i);
    else
        R_peaks(end+1) = PRT_peaks(i);
        R_locs(end+1) = PRT_locs(i);
    end
end

[QS_peaks,QS_locs] = findpeaks(-detrend_avg_healthy,MinPeakHeight=0.020,MinPeakProminence=0.03);
for i = 1:length(QS_peaks)
    if mod(i-2,3) == 0
        Q_peaks(end+1) = -QS_peaks(i);
        Q_locs(end+1) = QS_locs(i);
    elseif mod(i,3) == 0
        S_peaks(end+1) = -QS_peaks(i);
        S_locs(end+1) = QS_locs(i);
    end
end
    
figure(Name = 'PQRST Plot')
plot(healthy_time,detrend_avg_healthy,'-');
hold on
scatter(healthy_time(P_locs),P_peaks,'v','filled');
scatter(healthy_time(Q_locs),Q_peaks,'^','filled');
scatter(healthy_time(R_locs),R_peaks,'v','filled');
scatter(healthy_time(S_locs),S_peaks,'^','filled');
scatter(healthy_time(T_locs),T_peaks,'v','filled');
legend('','P','Q','R','S','T');
xlabel('Time (s)');
ylabel('mV');
title('Mean ECG Signal')

% Measured Heart Rate
QRS_int = [];
for i = 1:length(S_locs)-1
    QRS_int(end+1) = healthy_time(R_locs(i+1))-healthy_time(R_locs(i));
end
average_QRS_int = mean(QRS_int)
bpm = 60/average_QRS_int

% Maximum and Minimum
healthy_max = max(detrend_avg_healthy)
healthy_min = min(detrend_avg_healthy)

% Average Interval Calculations
average_PQ_int = mean(healthy_time(Q_locs)-healthy_time(P_locs))
average_PR_int = mean(healthy_time(R_locs)-healthy_time(P_locs))
average_QT_int = mean(healthy_time(T_locs(2:end))-healthy_time(Q_locs))

% MEA (Should probably check this part w/Jeffrey)
[peaks_healthy_I,~] = findpeaks(detrend_healthy_leads(:,1),MinPeakHeight=0.15,MinPeakDistance=20);
[peaks_healthy_III,~] = findpeaks(detrend_healthy_leads(:,3),MinPeakHeight=0.3,MinPeakDistance=20);
x1 = abs(mean(healthy_leadI))*cosd(0);
y1 = abs(mean(healthy_leadI))*sind(0);
x2 = abs(mean(healthy_leadIII))*cosd(120);
y2 = abs(mean(healthy_leadIII))*sind(120);
slope = tand(120);
slope_tang = -1/slope;
y3 = slope_tang*(x1-x2)+y2;
magnitude = sqrt(x1^2 + y3^2);
dir = atan2d(y3,x1);
% magnitude2 = sqrt((x1+x2)^2+(y1+y2)^2);
% dir2 = atan2d((y1+y2),(x1+x2));
figure(Name = 'Mean Axis of Depolarization')
c = compass([x1,x2,magnitude*cosd(dir)],[y1,y2,magnitude*sind(dir)]);
c(3).LineWidth = 2;
c(3).Color = 'r';
view(0,-90)
title('Mean Axis of Depolarization')
legend('Healthy Lead I','Healthy Lead III','MEA')
%% Part 3 (NOT done yet)
clear, clc, close all
diseased = readtable("Lab2_Disease_Data_ECG.xlsx");
diseased_time = diseased.Time;
a = diseased.LeadI;
b = diseased.LeadII;
c = diseased.LeadIII;
d = diseased.aVR;
e = diseased.aVL;
f = diseased.aVF;

diseased_leads = [a,b,c,d,e,f];
detrend_diseased_leads = [];
for i = 1:6
    [p, s, mu] = polyfit(diseased_time,diseased_leads(:,i),4);
    detrend_diseased_leads(:,end+1) = diseased_leads(:,i) - polyval(p, diseased_time, [], mu);
end

figure
subplot_titles_2 = ["Diseased Lead I","Diseased Lead II","Diseased Lead III","Diseased Lead aVR","Diseased Lead aVL","Diseased Lead aVF"];
for i = 1:6
    subplot(3,2,i)
    plot(diseased_time,detrend_diseased_leads(:,i))
    title(subplot_titles_2(i))
    xlabel('Time (s)')
    ylabel('mV')
end


% Average leads/detrend data
average_diseased_lead = mean(detrend_diseased_leads,2)
[p, s, mu] = polyfit(diseased_time,average_diseased_lead,4);
detrend_avg_diseased = average_diseased_lead - polyval(p, diseased_time, [], mu);

% Initialize PQRST arrays
P_peaks = [];
P_locs = [];
Q_peaks = [];
Q_locs = [];
T_peaks = [];
T_locs = [];
S_peaks = [];
S_locs = [];
R_peaks = [];
R_locs = [];

% Find PRT peaks
[PRT_peaks,PRT_locs] = findpeaks(detrend_avg_diseased,MinPeakHeight=0.01,MinPeakDistance=20);

for i = 1:length(PRT_peaks)
    if mod(i-1,3) == 0
        T_peaks(end+1) = PRT_peaks(i);
        T_locs(end+1) = PRT_locs(i);
    elseif  mod(i-2,3) == 0
        P_peaks(end+1) = PRT_peaks(i);
        P_locs(end+1) = PRT_locs(i);
    else
        R_peaks(end+1) = PRT_peaks(i);
        R_locs(end+1) = PRT_locs(i);
    end
end

[QS_peaks,QS_locs] = findpeaks(-detrend_avg_diseased,MinPeakHeight=0.020,MinPeakProminence=0.03);
for i = 1:length(QS_peaks)
    if mod(i-2,3) == 0
        Q_peaks(end+1) = -QS_peaks(i);
        Q_locs(end+1) = QS_locs(i);
    elseif mod(i,3) == 0
        S_peaks(end+1) = -QS_peaks(i);
        S_locs(end+1) = QS_locs(i);
    end
end
    
figure(Name = 'PQRST Plot')
p = plot(diseased_time,average_diseased_lead,'-');
hold on
scatter(diseased_time(P_locs),P_peaks,'v','filled');
scatter(diseased_time(Q_locs),Q_peaks,'^','filled');
scatter(diseased_time(R_locs),R_peaks,'v','filled');
scatter(diseased_time(S_locs),S_peaks,'^','filled');
scatter(diseased_time(T_locs),T_peaks,'v','filled');
legend('','P','Q','R','S','T');
xlabel('Time (s)');
ylabel('mV');

% Measured Heart Rate
QRS_int = [];
for i = 1:length(R_locs)-1
    QRS_int(end+1) = diseased_time(R_locs(i+1))-diseased_time(R_locs(i));
end
average_QRS_int = mean(QRS_int)
bpm = 60/average_QRS_int

% Maximum and Minimum
healthy_max = max(detrend_avg_diseased)
healthy_min = min(detrend_avg_diseased)

% Average Interval Calculations
average_PQ_int = mean(diseased_time(Q_locs)-diseased_time(P_locs))
average_PR_int = mean(diseased_time(R_locs)-diseased_time(P_locs))
average_QT_int = mean(diseased_time(T_locs(2:end))-diseased_time(Q_locs))