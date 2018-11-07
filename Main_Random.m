clear%vars -except TarTruePathStruct3 TarTruePathStruct4;
close all;
delete AllData1.mat
delete AllData2.mat
delete AllData3.mat
delete AllData4.mat
load('MaxAcc.mat');
for caseLoop  = ["C","F","I","II"]
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
            UKF_TP = zeros(length(clutterVersion),nYaw-nYawStart);
            UKF_RMS = zeros(length(clutterVersion),nYaw-nYawStart);
            save('performance','UKF_TP','UKF_RMS');        
        end
        for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
            totalSim = 2000;        
            temp.UKFTP = 0;
            temp.UKFRMS = 0;

            for nSim = 1:totalSim
                if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                    save('loopParam','clutterVersion','totalSim')                   
                else                
                    load('loopParam');
                end
                if nSim == totalSim
                    load('performance');
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
                transWave.true = Planner_True(ParPath,Timings.impInstVec);
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                acc = MaxAcc(2,(MaxAcc(1,:) == yawDevRateLoop));
                noiseAcc = [acc;acc];
                [FilterEstimateUKF,transWave.estUKF,Measurements] = UKF_PDAF(ParPath,ParGen,Timings,Measurements,noiseAcc,ClutterStruct,caseLoop,ParStdDev);

                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                impLength = length(Timings.impInstVec);
                Threshold.cw = 250;
                Threshold.fm = 250;                
                
                XMinusVecUKF = FilterEstimateUKF.XkMinusVec';
                diffBtwTrueExpUKF = (ParPath.TarPathMat(Timings.impInstVec,:) - XMinusVecUKF(:,1:2)).^2;
                diffBtwTrueExpUKF = sqrt(diffBtwTrueExpUKF(:,1) + diffBtwTrueExpUKF(:,2)); 
                CwIndx = find(transWave.estUKF == "C");
                FmIndx = find(transWave.estUKF == "F");
                trackPointsUKF = sum(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);
                temp.UKFTP = temp.UKFTP + (trackPointsUKF/impLength)*100;
                trackCommInd = [FmIndx(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw)];
                temp.UKFRMS = temp.UKFRMS + sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)/length(trackCommInd));
                if nSim == totalSim                
                    UKF_TP(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = temp.UKFTP/totalSim;
                    UKF_RMS(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = temp.UKFRMS/totalSim;
                end
   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                save('loopParam','clutterVersion','totalSim');
                if nSim == totalSim
                    save('performance','UKF_TP','UKF_RMS');
                    clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                        yawDevRateLoopJump nYawStart caseLoop  MaxAcc;
                else
                    clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                        yawDevRateLoopJump nYawStart temp caseLoop MaxAcc;
                end

            end
            disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
        end
        disp(['cluttterLoop:',num2str(clutterLoop)]);
    end
    switch caseLoop
        case "C"
            load('performance');            
            save('AllData1','UKF_TP','UKF_RMS','caseLoop');
        case "F"
            load('performance');
            save('AllData2','UKF_TP','UKF_RMS','caseLoop');
        case "I"
            load('performance');
            save('AllData3','UKF_TP','UKF_RMS','caseLoop');
        case "II"
            load('performance');
            save('AllData4','UKF_TP','UKF_RMS','caseLoop');
    end
    disp(["caseLoop:" caseLoop]);
    clearvars -except MaxAcc ;
end 
