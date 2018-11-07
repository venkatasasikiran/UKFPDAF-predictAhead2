close all;
clear;
load("MaxAcc.mat")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
turnRate = 5;
ParPath = Struct_Param_Path_Init(true,turnRate,90,[6000,4000],400,100);
% ParPath.nWaves = 1;
% ParPath.maxDepth = 2;
ParGen = Struct_Param_Gen_Init(1/100000);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ParPath,~]= Path_Target_With_Models(ParPath,ParGen);
figure(1);
plot(ParPath.TarPathMat(:,1),ParPath.TarPathMat(:,2));
title('Ground true');
hold off; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Timings] = Timing_Model(ParPath,ParGen);
% if length(Timings.impInstVec)< 30
%    while(length(Timings.impInstVec)< 30)
%        [ParPath,devTimeVec,~]= Path_Target(ParPath,ParGen);
%        [Timings] = Timing_Model(ParPath,ParGen);
%    end
% end
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParStdDev = Struct_Param_Std_Dev([40,10],deg2rad(1),.5);
[Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ClutterStruct = Clutter_Generation(ParPath.TarMaxRho+150,ParGen.clutterDensity,Timings.impInstVec);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transWave.true = Planner_True(ParPath,Timings.impInstVec);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc = MaxAcc(2,(MaxAcc(1,:) == turnRate));
noiseAcc = [acc;acc];
% [FilterEstimateUKF,transWave.estUKF,Measurements] = UKF_PDAF(ParPath,ParGen,Timings,Measurements,...
%                                                                 noiseAcc,ClutterStruct,"I",ParStdDev);
[FilterEstimateUKF,transWave.estUKF,Measurements] =UKF_PDAF3(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,ParStdDev);                                                            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
impLength = length(Timings.impInstVec);
Threshold.cw = 250;
Threshold.fm = 250;                
                
% XMinusVecUKF = FilterEstimateUKF.XkMinusVec';
XMinusVecUKF = FilterEstimateUKF.XkPlusVec';

diffBtwTrueExpUKF = (ParPath.TarPathMat(Timings.impInstVec(1:end-2),:) - XMinusVecUKF(2:end-2,1:2)).^2;
diffBtwTrueExpUKF = sqrt(diffBtwTrueExpUKF(:,1) + diffBtwTrueExpUKF(:,2)); 
CwIndx = find(transWave.estUKF(1:end-2) == "C");
FmIndx = find(transWave.estUKF(1:end-2) == "F");
trackPointsUKF = sum(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);
UKFTP = (trackPointsUKF/(impLength-2))*100;

trackCommInd = [FmIndx(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw)];
UKFRMS       = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)/length(trackCommInd));
  
display(UKFTP);
display(UKFRMS);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure(2)
plot(FilterEstimateUKF.XkMinusVec(1,:),FilterEstimateUKF.XkMinusVec(2,:),'r-*');
hold on;        
plot(ParPath.TarPathMat(Timings.impInstVec,1),ParPath.TarPathMat(Timings.impInstVec,2),'-o');
title('Filtered and True UKF');
hold off;  
            
figure(3);
plot(diffBtwTrueExpUKF,'-*');
title('diff true exp');

figure(4);
EstPosCwIndx = find(transWave.estUKF == "C");
EstPosFmIndx = find(transWave.estUKF == "F");
plot(FilterEstimateUKF.XkMinusVec(1,:),FilterEstimateUKF.XkMinusVec(2,:));
hold on;
plot(FilterEstimateUKF.XkMinusVec(1,EstPosCwIndx),FilterEstimateUKF.XkMinusVec(2,EstPosCwIndx),'ro');
plot(FilterEstimateUKF.XkMinusVec(1,EstPosFmIndx),FilterEstimateUKF.XkMinusVec(2,EstPosFmIndx),'bd');
hold off;
title("red-circle: CW blue-diamond: FM")

figure(5)  
subplot(2,2,1)
plot(Timings.impInstVec,ParPath.TrueVelocities(:,1));
title('true x vel')

subplot(2,2,3)
plot(Timings.impInstVec,ParPath.TrueVelocities(:,2));
title('true y vel')

subplot(2,2,2)
plot(Timings.impInstVec,FilterEstimateUKF.XkMinusVec(3,:));
title('Estimated x vel UKF')

subplot(2,2,4)
plot(Timings.impInstVec,FilterEstimateUKF.XkMinusVec(4,:));
title('Estimated y vel UKF')  


figure(6)
polarplot(ParPath.TarPathMatInPol(:,1),ParPath.TarPathMatInPol(:,2));
hold on;
[thetaEst,rhoEst] = cart2pol(FilterEstimateUKF.XkMinusVec(1,:),FilterEstimateUKF.XkMinusVec(2,:));
polarplot(thetaEst,rhoEst);
polarscatter(thetaEst(EstPosCwIndx),rhoEst(EstPosCwIndx),'red','o','filled');
polarscatter(thetaEst(EstPosFmIndx),rhoEst(EstPosFmIndx),'blue','d','filled');
hold off;
title("red-circle: CW blue-diamond: FM")

