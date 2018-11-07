clearvars -except MeasStruct2 TarTruePathStruct2;
close all;
% load('PathSet2.mat');
load('MaxAcc.mat');
% load('MeasStructs.mat')


for caseLoop  = ["C","F","I","II"]
    clutterVersion = 100000;
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
            totalSim = 100;        
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
                if yawDevRateLoop <= 2.5
                    structX = (yawDevRateLoop/yawDevRateLoopJump)+1 - yawDevRateLoopMin/yawDevRateLoopJump;
                    ParPath.TarPathMat = TarTruePathStruct1(structX,nSim).TarPathMat;
                    ParPath.TarPathMatInPol = TarTruePathStruct1(structX,nSim).TarPathMatInPol;
                    ParPath.TarMaxRho = TarTruePathStruct1(structX,nSim).TarMaxRho;
                else
                    if yawDevRateLoop > 4.5
                        error('change the code and choose the appropriate pathset')
                    else
                        structX = yawDevRateLoop/yawDevRateLoopJump - 3/yawDevRateLoopJump + 1;
                        ParPath.TarPathMat = TarTruePathStruct2(structX,nSim).TarPathMat;
                        ParPath.TarPathMatInPol = TarTruePathStruct2(structX,nSim).TarPathMatInPol;
                        ParPath.TarMaxRho = TarTruePathStruct2(structX,nSim).TarMaxRho;
                    end
                end

%              
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Timings] = Timing_Model(ParPath,ParGen);
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).impInstVec = Timings.impInstVec;          
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ParStdDev = Struct_Param_Std_Dev([40,10],deg2rad(1),.5);              
                if yawDevRateLoop <= 2.5
                    Measurements = MeasStruct1(structX,nSim).Measurements;
                    ParPath.TrueVelocities = MeasStruct1(structX,nSim).TrueVel ;
                else
                    Measurements = MeasStruct2(structX,nSim).Measurements;
                    ParPath.TrueVelocities = MeasStruct2(structX,nSim).TrueVel ;
                end
             

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
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).XMinusVecUKF = XMinusVecUKF;

                diffBtwTrueExpUKF = (ParPath.TarPathMat(Timings.impInstVec,:) - XMinusVecUKF(:,1:2)).^2;
                diffBtwTrueExpUKF = sqrt(diffBtwTrueExpUKF(:,1) + diffBtwTrueExpUKF(:,2)); 
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).diffBtwTrueExpUKF = diffBtwTrueExpUKF;
                CwIndx = find(transWave.estUKF == "C");
                FmIndx = find(transWave.estUKF == "F");
                trackPointsUKF = sum(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).UKFTP = (trackPointsUKF/impLength)*100;
                temp.UKFTP = temp.UKFTP + (trackPointsUKF/impLength)*100;

                trackCommInd = [FmIndx(diffBtwTrueExpUKF(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExpUKF(CwIndx) <= Threshold.cw)];
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).trackCommInd = trackCommInd;


                temp.UKFRMS = temp.UKFRMS + sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)/length(trackCommInd));
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).UKFRMS = sqrt(sum(diffBtwTrueExpUKF(trackCommInd).^2)/length(trackCommInd));
                if nSim == totalSim                
                    UKF_TP(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = temp.UKFTP/totalSim;
                    UKF_RMS(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = temp.UKFRMS/totalSim;
                end
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).turnRate = yawDevRateLoop;
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).nSim = nSim;
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).clutterDensity = 1/clutterVersion(clutterLoop);
                AllData(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).caseLoop = caseLoop;        

   


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                save('loopParam','clutterVersion','totalSim');
                if nSim == totalSim
                    save('performance','UKF_TP','UKF_RMS');
                    clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                        yawDevRateLoopJump nYawStart caseLoop  MaxAcc TarTruePathStruct1 ... 
                        TarTruePathStruct2 AllData MeasStruct1 MeasStruct2;
                else
                    clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                        yawDevRateLoopJump nYawStart temp caseLoop MaxAcc TarTruePathStruct1 ...
                        TarTruePathStruct2 AllData  MeasStruct1 MeasStruct2;
                end

            end
            disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
        end
        disp(['cluttterLoop:',num2str(clutterLoop)]);
    end
    switch caseLoop
        case "C"
            load('performance');
            AllData1 = AllData;
            save('AllData1','AllData1','-v7.3','UKF_TP','UKF_RMS','caseLoop');
%             save('AllData1','UKF_TP','UKF_RMS','caseLoop');
        case "F"
            load('performance');
            AllData2 = AllData;
            save('AllData2','AllData2','-v7.3','UKF_TP','UKF_RMS','caseLoop');
%             save('AllData2','UKF_TP','UKF_RMS','caseLoop');
        case "I"
            load('performance');
            AllData3 = AllData;
            save('AllData3','AllData3','-v7.3','UKF_TP','UKF_RMS','caseLoop');
%             save('AllData3','UKF_TP','UKF_RMS','caseLoop');
        case "II"
            load('performance');
            AllData4 = AllData;
            save('AllData4','AllData4','-v7.3','UKF_TP','UKF_RMS','caseLoop');
    end
    disp(["caseLoop:" caseLoop]);
    clearvars -except MaxAcc TarTruePathStruct1 TarTruePathStruct2 MeasStruct1 MeasStruct2;

end