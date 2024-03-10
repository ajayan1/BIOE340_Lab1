%% Part 1 (done)
clear; clc; close all

% Import healthy data
healthy = readtable("Lab2_Healthy_Data_ECG.xlsx");

healthy_time = healthy.Time;
healthy_leadI = healthy.LeadI./1000;
healthy_leadII = healthy.LeadII./1000;
healthy_leadIII = healthy.LeadIII./1000;
healthy_aVR = healthy.aVR./1000;
healthy_aVL = healthy.aVL./1000;
healthy_aVF = healthy.aVF./1000;

% Detrend healthy data
healthy_leads = [healthy_leadI,healthy_leadII,healthy_leadIII,healthy_aVR,healthy_aVL,healthy_aVF];
detrend_healthy_leads = [];
for i = 1:6
    [p, ~, mu] = polyfit(healthy_time,healthy_leads(:,i),7);
    detrend_healthy_leads(:,end+1) = healthy_leads(:,i) - polyval(p, healthy_time, [], mu);
end

% Plot ECG data
figure(Name = '6-Lead ECG Healthy')
subplot_titles_1 = ["Healthy Lead I","Healthy Lead II","Healthy Lead III","Healthy Lead aVR","Healthy Lead aVL","Healthy Lead aVF"];
for i = 1:6
    subplot(2,3,i)
    plot(healthy_time,detrend_healthy_leads(:,i))
    title(subplot_titles_1(i))
    xlabel('Time (s)')
    ylabel('mV')
end

% Derive LeadIII, aVF, aVL, aVR
Cal_LeadIII = detrend_healthy_leads(:,2) - detrend_healthy_leads(:,1);
Cal_aVF = ((2.*detrend_healthy_leads(:,2))-detrend_healthy_leads(:,1))./(sqrt(3));
Cal_aVL = ((2.*detrend_healthy_leads(:,1))-detrend_healthy_leads(:,2))./(sqrt(3));
Cal_aVR = -(detrend_healthy_leads(:,2)+detrend_healthy_leads(:,1))./(sqrt(3));

a = isequal(Cal_LeadIII,detrend_healthy_leads(:,3))

% Compare derived leads to actual
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
%% Part 2 (ask jeffrey about MEA)

% Average leads/detrend data
average_healthy_lead = mean(detrend_healthy_leads,2);
[p, ~, mu] = polyfit(healthy_time,average_healthy_lead,7);
detrend_avg_healthy = average_healthy_lead - polyval(p, healthy_time, [], mu);
smoothECG_healthy = sgolayfilt(detrend_avg_healthy,7,21);

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

% Find P, R, and T
[PRT_peaks,PRT_locs] = findpeaks(smoothECG_healthy,NPeaks=25,MinPeakHeight=0.01,MinPeakDistance=20);

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

% Find Q and S
[QS_peaks,QS_locs] = findpeaks(-smoothECG_healthy,MinPeakHeight=0.020,MinPeakProminence=0.03);
for i = 1:length(QS_peaks)
    if mod(i-2,3) == 0
        Q_peaks(end+1) = -QS_peaks(i);
        Q_locs(end+1) = QS_locs(i);
    elseif mod(i,3) == 0
        S_peaks(end+1) = -QS_peaks(i);
        S_locs(end+1) = QS_locs(i);
    end
end
    
% Plot PQRTS
figure(Name = 'PQRST Plot')
plot(healthy_time,smoothECG_healthy,'-');
hold on
scatter(healthy_time(P_locs),P_peaks,'v','filled');
scatter(healthy_time(Q_locs),Q_peaks,'^','filled');
scatter(healthy_time(R_locs),R_peaks,'v','filled');
scatter(healthy_time(S_locs),S_peaks,'^','filled');
scatter(healthy_time(T_locs),T_peaks,'v','filled');
legend('','P','Q','R','S','T');
xlabel('Time (s)');
ylabel('mV');
title('Mean ECG Signal Healthy')

% Measure Heart Rate
QRS_int = [];
for i = 1:length(S_locs)-1
    QRS_int(end+1) = healthy_time(R_locs(i+1))-healthy_time(R_locs(i));
end
average_QRS_int = mean(QRS_int)
bpm = 60/average_QRS_int

% Maximum and Minimum
healthy_max = max(smoothECG_healthy)
healthy_min = min(smoothECG_healthy)

% Average Interval Calculations
average_PQ_int = mean(healthy_time(Q_locs)-healthy_time(P_locs))
average_PR_int = mean(healthy_time(R_locs)-healthy_time(P_locs))
average_QT_int = mean(healthy_time(T_locs(2:end))-healthy_time(Q_locs))

% MEA (Should probably check this part w/Jeffrey)
[peaks_healthy_I,~] = findpeaks(detrend_healthy_leads(:,1),MinPeakHeight=0.15,MinPeakDistance=20);
[peaks_healthy_III,~] = findpeaks(detrend_healthy_leads(:,3),MinPeakHeight=0.3,MinPeakDistance=20);
x1 = abs(mean(peaks_healthy_I))*cosd(0);
y1 = abs(mean(peaks_healthy_I))*sind(0);
x2 = abs(mean(peaks_healthy_III))*cosd(120);
y2 = abs(mean(peaks_healthy_III))*sind(120);
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
title('Mean Axis of Depolarization Healthy')
legend('Healthy Lead I','Healthy Lead III','MEA')
%% Part 3 (need to ask jeffrey about)

diseased = readtable("Lab2_Disease_Data_ECG.xlsx");
diseased_time = diseased.Time;
diseased_LeadI = diseased.LeadI;
diseased_LeadII = diseased.LeadII;
diseased_LeadIII = diseased.LeadIII;
diseased_aVR = diseased.aVR;
diseased_aVL = diseased.aVL;
diseased_aVF = diseased.aVF;

diseased_leads = [diseased_LeadI,diseased_LeadII,diseased_LeadIII,diseased_aVR,diseased_aVL,diseased_aVF];
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
average_diseased_lead = mean(detrend_diseased_leads,2);
[p, s, mu] = polyfit(diseased_time,average_diseased_lead,7);
detrend_avg_diseased = average_diseased_lead - polyval(p, diseased_time, [], mu);
smoothECG_diseased = sgolayfilt(detrend_avg_diseased,7,21);

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
[PRT_peaks,PRT_locs] = findpeaks(smoothECG_diseased,MinPeakHeight=0.01,MinPeakDistance=20);

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

[QS_peaks,QS_locs] = findpeaks(-smoothECG_diseased,MinPeakHeight=0.020,MinPeakProminence=0.03);
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
plot(diseased_time,smoothECG_diseased,'-');
hold on
scatter(diseased_time(P_locs),P_peaks,'v','filled');
scatter(diseased_time(Q_locs),Q_peaks,'^','filled');
scatter(diseased_time(R_locs),R_peaks,'v','filled');
scatter(diseased_time(S_locs),S_peaks,'^','filled');
scatter(diseased_time(T_locs),T_peaks,'v','filled');
legend('','P','Q','R','S','T');
xlabel('Time (s)');
ylabel('mV');
xlabel('Time (s)');
ylabel('mV');
title('Mean ECG Signal Diseased')

% Measured Heart Rate
QRS_int = [];
for i = 1:length(R_locs)-1
    QRS_int(end+1) = diseased_time(R_locs(i+1))-diseased_time(R_locs(i));
end
average_QRS_int = mean(QRS_int)
bpm = 60/average_QRS_int

% Maximum and Minimum
diseased_max = max(smoothECG_diseased)
diseased_min = min(smoothECG_diseased)

% Average Interval Calculations (Probably cant do this)
average_PQ_int = mean(diseased_time(Q_locs)-diseased_time(P_locs))
average_PR_int = mean(diseased_time(R_locs)-diseased_time(P_locs))
average_QT_int = mean(diseased_time(T_locs(2:end))-diseased_time(Q_locs))

% MEA (Should probably check this part w/Jeffrey)
[peaks_diseased_I,~] = findpeaks(detrend_diseased_leads(:,1),MinPeakHeight=0.15,MinPeakDistance=20);
[peaks_diseased_III,~] = findpeaks(detrend_diseased_leads(:,3),MinPeakHeight=0.3,MinPeakDistance=20);
x1 = abs(mean(peaks_diseased_I))*cosd(0);
y1 = abs(mean(peaks_diseased_I))*sind(0);
x2 = abs(mean(peaks_diseased_III))*cosd(120);
y2 = abs(mean(peaks_diseased_III))*sind(120);
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
title('Mean Axis of Depolarization Diseased')
legend('Healthy Lead I','Healthy Lead III','MEA')

%% Diagnosis (maybe ok)
[peaks1,locs1] = findpeaks(detrend_diseased_leads(:,1),MinPeakHeight=0.15,MinPeakDistance=20);
[peaks2,locs2] = findpeaks(detrend_diseased_leads(:,2),MinPeakHeight=0.15,MinPeakDistance=20);
[peaks3,locs3] = findpeaks(detrend_diseased_leads(:,3),MinPeakHeight=0.15,MinPeakDistance=20);
[npeaks1,nlocs1] = findpeaks(-detrend_diseased_leads(:,1),MinPeakHeight=0.15,MinPeakDistance=20);
[npeaks2,nlocs2] = findpeaks(-detrend_diseased_leads(:,2),MinPeakHeight=0.15,MinPeakDistance=20);
[npeaks3,nlocs3] = findpeaks(-detrend_diseased_leads(:,3),MinPeakHeight=0.15,MinPeakDistance=20);

v1 = [];
v2 = [];
v3 = [];
for i = 1:length(peaks1)
    if mod(i,3) == 0
        v1(end+1) = peaks1(i) - npeaks1(i);
        v2(end+1) = peaks2(i) - npeaks2(i);
        v3(end+1) = peaks3(i) - npeaks3(i);
    end
end

sum_QRS_voltage = mean(v1) + mean(v2) + mean(v3)
average_QRS_int = mean(diseased_time(S_locs)-diseased_time(Q_locs))

issues = string();
possible_diseases = string();

if sum_QRS_voltage < 0.5
    issues(end+1) = "Low Voltage";
    possible_diseases(end+1) = "Pericardial fluid buildup";
    possible_diseases(end+1) = "Pulmonary emphysema";
    possible_diseases(end+1) = "Previous myocardial infarctions/diminished cardiac muscle mass";
elseif sum_QRS_voltage >2.0
    issues(end+1) = "High Voltage";
    possible_diseases(end+1) = "Hypertrophy (High Voltage)";
end
if average_QRS_int > 0.08
    issues(end+1) = "Prolonged QRS Wave";
    if average_QRS_int <= 0.12
        possible_diseases(end+1) = "Hypertrophy (Prolonged QRS Wave)";
        possible_diseases(end+1) = "Dilation";
    elseif average_QRS_int > 0.12
        possible_diseases(end+1) = "Damage to cardiac muscle";
        possible_diseases(end+1) = "Blocks in the Purkinje system";
    end
end

if sum_QRS_voltage >= 0.5 && sum_QRS_voltage <= 2.0 && average_QRS_int >= 0.06 && average_QRS_int <= 0.08
    fprintf('Patient is healthy.\n')
else
    fprintf('Issues:\n')
    fprintf('%s\n',issues(2:end))
    fprintf('\nPossible Diseases:\n')
    fprintf('%s\n',possible_diseases(2:end))
end
