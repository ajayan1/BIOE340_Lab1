Healty = readtable("Lab2_Healty_Data_ECG.xlsx");
Diseased = readtable("Lab2_Disease_Data_ECG.xlsx");

Healty_Time = Healty.Time;
Healty_LeadI = Healty.LeadI;
Healty_LeadII = Healty.LeadII;
Healty_LeadIII = Healty.LeadIII;
Healty_aVR = Healty.aVR;
Healty_aVL = Healty.aVL;
Healty_aVF = Healty.aVF;

Diseased_Time = Diseased.Time;
Diseased_LeadI = Diseased.LeadI;
Diseased_LeadII = Diseased.LeadII;
Diseased_LeadIII = Diseased.LeadIII;
Diseased_aVR = Diseased.aVR;
Diseased_aVL = Diseased.aVL;
Diseased_aVF = Diseased.aVF;

figure(1)
subplot(6,1,1)
plot(Healty_Time,Healty_LeadI)
title('Healthy Lead I')

subplot(6,1,2)
plot(Healty_Time,Healty_LeadII)
title('Healthy Lead II')

subplot(6,1,3)
plot(Healty_Time,Healty_LeadIII)
title('Healthy Lead III')

subplot(6,1,4)
plot(Healty_Time,Healty_aVR)
title('Healthy Lead aVR')

subplot(6,1,5)
plot(Healty_Time,Healty_aVL)
title('Healthy Lead aVL')

subplot(6,1,6)
plot(Healty_Time,Healty_aVF)
title('Healthy Lead aVF')


Cal_LeadIII = Healty_LeadII - Healty_LeadI;
Cal_aVF = ((2.*Healty_LeadII)-Healty_LeadI)./(sqrt(3));
Cal_aVL = ((2.*Healty_LeadI)-Healty_LeadII)./(sqrt(3));
Cal_aVR = -(Healty_LeadII+Healty_LeadI)./(sqrt(3));


a = isequal(Cal_LeadIII,Healty_LeadIII)

figure(2)

subplot(4,2,1)
plot(Healty_Time,Cal_LeadIII)
title('Derived Lead III')

subplot(4,2,2)
plot(Healty_Time,Healty_LeadIII)
title('Healthy Lead III')

subplot(4,2,3)
plot(Healty_Time,Cal_aVR)
title('Derived Lead aVR')

subplot(4,2,4)
plot(Healty_Time,Healty_aVR)
title('Healthy Lead aVR')

subplot(4,2,5)
plot(Healty_Time,Cal_aVF)
title('Derived Lead aVF')

subplot(4,2,6)
plot(Healty_Time,Healty_aVF)
title('Healthy Lead aVF')

subplot(4,2,7)
plot(Healty_Time,Cal_aVL)
title('Derived Lead aVL')

subplot(4,2,8)
plot(Healty_Time,Healty_aVL)
title('Healthy Lead aVL')


