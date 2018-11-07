clear;
close all;
load('MaxAcc.mat');
clutterVersion = [100000];
for clutterLoop = 1:length(clutterVersion)
    yawDevRateLoopMin = 4;
    yawDevRateLoopJump = .5;
    yawDevRateLoopMax = 4.5;
    nYaw = (yawDevRateLoopMax/yawDevRateLoopJump)+1;
    nYawStart = yawDevRateLoopMin/yawDevRateLoopJump;
    if clutterLoop == 1
        delete performance.mat
        delete loopParam.mat                   
        UKF_TP_C = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS_C = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_TP_F = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS_F = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_TP_I = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS_I = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_TP_II = zeros(length(clutterVersion),nYaw-nYawStart);
        UKF_RMS_II = zeros(length(clutterVersion),nYaw-nYawStart);
%         tempValues(nYaw-nYawStart) = struct;
        save('performance','UKF_TP_C','UKF_RMS_C','UKF_TP_F','UKF_RMS_F',...
             'UKF_TP_I','UKF_RMS_I','UKF_TP_II','UKF_RMS_II');        
    end
    for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
    totalSim = 2000;        
    temp.UKFTP_C = zeros(totalSim,1);
    temp.UKFRMS_C = zeros(totalSim,1);
    temp.UKFTP_F = zeros(totalSim,1);
    temp.UKFRMS_F = zeros(totalSim,1);
    temp.UKFTP_I = zeros(totalSim,1);
    temp.UKFRMS_I = zeros(totalSim,1);
    temp.UKFTP_II = zeros(totalSim,1);
    temp.UKFRMS_II = zeros(totalSim,1);
        for nSim = 1:totalSim
            if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                save('loopParam','clutterVersion','totalSim')                   
            else                
                load('loopParam');
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],400,100);
            ParPath.maxDepth = 2;
            ParGen = Struct_Param_Gen_Init(1/clutterVersion(clutterLoop));
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
            [ParPath,~]= Path_Target_With_Models(ParPath,ParGen);
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Timings] = Timing_Model(ParPath,ParGen);
            if length(Timings.impInstVec)< 30
               while(length(Timings.impInstVec)< 30)
                   [ParPath,~]= Path_Target_With_Models(ParPath,ParGen);
                   [Timings] = Timing_Model(ParPath,ParGen);
               end
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ParStdDev = Struct_Param_Std_Dev([40,10],deg2rad(1),.5);                  
            [Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev); 
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ClutterStruct = Clutter_Generation(ParPath.TarMaxRho,ParGen.clutterDensity,Timings.impInstVec);
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 transWave.true = Planner_True(ParPath,Timings.impInstVec);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            acc = MaxAcc(2,(MaxAcc(1,:) == yawDevRateLoop));
            noiseAcc = [acc;acc];
            impLength = length(Timings.impInstVec);
            Threshold.cw = 250;
            Threshold.fm = 250;
            for caseLoop = ["C","F","I","II"]
                [FilterEstimateUKF,transWave.estUKF,Measurements] = UKF_PDAF(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,caseLoop,ParStdDev);

                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           

                XMinusVecUKF = FilterEstimateUKF.XkMinusVec';
                diffBtwTrueExpUKF = (ParPath.TarPathMat(Timings.impInstVec,:) - XMinusVecUKF(:,1:2)).^2;
                diffBtwTrueExpUKF = sqrt(diffBtwTrueExpUKF(:,1) + diffBtwTrueExpUKF(:,2)); 
                CwIndx = find(transWave.estUKF == "C");
                FmIndx = find(transWave.estUKF == "F");
                trackPointsUKF = sum(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);
                trackCommInd = [FmIndx(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw)];
                switch caseLoop
                    case "C" 
                        temp.UKFTP_C(nSim) = (trackPointsUKF/impLength)*100;        
                        temp.UKFRMS_C(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "F"
                        temp.UKFTP_F(nSim) = (trackPointsUKF/impLength)*100;        
                        temp.UKFRMS_F(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "I"
                        temp.UKFTP_I(nSim) = (trackPointsUKF/impLength)*100;        
                        temp.UKFRMS_I(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "II"
                        temp.UKFTP_II(nSim) = (trackPointsUKF/impLength)*100;        
                        temp.UKFRMS_II(nSim) = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)...
                                                /length(trackCommInd));
                end
                
                if nSim == totalSim
                    load('performance');
                    UKF_TP_C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFTP_C);
                    UKF_RMS_C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFRMS_C);
                    UKF_TP_F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFTP_F);
                    UKF_RMS_F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFRMS_F);
                    UKF_TP_I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFTP_I);
                    UKF_RMS_I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFRMS_I);
                    UKF_TP_II(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFTP_II);
                    UKF_RMS_II(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(temp.UKFRMS_II);
                    tempValues(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = temp;
                end
            end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            save('loopParam','clutterVersion','totalSim');
            if nSim == totalSim
                save('performance','UKF_TP_C','UKF_RMS_C','UKF_TP_F','UKF_RMS_F',...
             'UKF_TP_I','UKF_RMS_I','UKF_TP_II','UKF_RMS_II','tempValues');
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                    yawDevRateLoopJump nYawStart  MaxAcc;
            else
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                    yawDevRateLoopJump nYawStart temp  MaxAcc;
            end

        end
        disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
    end
    disp(['cluttterLoop:',num2str(clutterLoop)]);    
    clearvars -except MaxAcc ;
end

