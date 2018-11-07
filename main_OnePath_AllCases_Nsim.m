close all;
clear;
load('MaxAcc.mat');
clutterVersion = [100000];
for clutterLoop = 1:length(clutterVersion)
    yawDevRateLoopMin = 4.5;
    yawDevRateLoopJump = .5;
    yawDevRateLoopMax = 7;
    nYaw = (yawDevRateLoopMax/yawDevRateLoopJump)+1;
    nYawStart = yawDevRateLoopMin/yawDevRateLoopJump;
    if clutterLoop == 1
        delete performance.mat
        delete loopParam.mat                   
       UKF_TP.C = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.C = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_TP.F = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.F = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_TP.I = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.I = zeros(length(clutterVersion),nYaw-nYawStart);       
       UKF_TP.IE = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.IE = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_TP.P = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.P = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_TP.PE = zeros(length(clutterVersion),nYaw-nYawStart);
       UKF_RMS.PE = zeros(length(clutterVersion),nYaw-nYawStart);
       save('performance','UKF_TP','UKF_RMS');        
    end
    for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
        totalSim = 1000;        
        tempTP.C = zeros(totalSim,1);
        tempRMS.C = zeros(totalSim,1);    
        tempTP.F = zeros(totalSim,1); 
        tempRMS.F = zeros(totalSim,1);
        tempTP.I = zeros(totalSim,1);
        tempRMS.I = zeros(totalSim,1);        
        tempTP.IE = zeros(totalSim,1);
        tempRMS.IE = zeros(totalSim,1);
        tempTP.P = zeros(totalSim,1);
        tempRMS.P = zeros(totalSim,1);
        tempTP.PE = zeros(totalSim,1);
        tempRMS.PE = zeros(totalSim,1);
        for nSim = 1:totalSim
            if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                save('loopParam','clutterVersion','totalSim')                   
            else                
                load('loopParam');                
            end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
%             ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],400,100);
            ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],400,100);
            ParGen = Struct_Param_Gen_Init(1/clutterVersion(clutterLoop));% clutterdensity inside the parenthesis

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ParPath,~]= Path_Target_With_Models(ParPath,ParGen);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Timings] = Timing_Model(ParPath,ParGen);                               
    %% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ParStdDev = Struct_Param_Std_Dev([40,10],deg2rad(1),.5);
            [Measurements,ParPath] = Measurement_Model(ParPath,Timings.impInstVec,ParStdDev);                                                            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ClutterStruct = Clutter_Generation(ParPath.TarMaxRho+150,ParGen.clutterDensity,Timings.impInstVec);
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            transWave.true = Planner_True(ParPath,Timings.impInstVec);
            impLength = length(Timings.impInstVec);
            Threshold.cw = 250;
            Threshold.fm = 250;
            
            acc = MaxAcc(2,(MaxAcc(1,:) == yawDevRateLoop));
            noiseAcc = [acc;acc];
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            for caseLoop = ["C","F","I","P","IE","PE"]
                if caseLoop == "P"
                    [FilterEstimate,transWave.estUKF,Measurements] =UKF_PDAF3(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,ParStdDev);
                    
                
                elseif caseLoop == "PE"
                    [FilterEstimate,transWave.estUKF,Measurements] =UKF_PDAF_PE(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,ParStdDev);
                else
                    [FilterEstimate,transWave.estUKF,Measurements] =UKF_PDAF(ParPath,ParGen,Timings,Measurements,...
                                                                        noiseAcc,ClutterStruct,caseLoop,ParStdDev);
                    
                end

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                XPlusVecUKF = FilterEstimate.XkPlusVec(:,2:end)';
                diffBtwTrueExp = (ParPath.TarPathMat(Timings.impInstVec,:) - XPlusVecUKF(:,1:2)).^2;
                diffBtwTrueExp = sqrt(diffBtwTrueExp(:,1) + diffBtwTrueExp(:,2)); 
                CwIndx = find(transWave.estUKF == "C");
                FmIndx = find(transWave.estUKF == "F");
                trackPoints = sum(diffBtwTrueExp(CwIndx) <= Threshold.cw) + sum(diffBtwTrueExp(FmIndx) <= Threshold.fm);
                trackCommInd = [FmIndx(diffBtwTrueExp(FmIndx) <= Threshold.fm);CwIndx(diffBtwTrueExp(CwIndx) <= Threshold.cw)];    
                switch caseLoop
                    case "C" 
                        tempTP.C(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.C(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "F"
                        tempTP.F(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.F(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "I"
                        tempTP.I(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.I(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "P"
                        tempTP.P(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.P(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "IE"                                            
                        tempTP.IE(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.IE(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                    case "PE"
                        tempTP.PE(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.PE(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                end
                if nSim == totalSim
                    load('performance');
                   UKF_TP.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.C);
                   UKF_RMS.C(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.C);
                   UKF_TP.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.F);
                   UKF_RMS.F(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.F);
                   UKF_TP.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.I);
                   UKF_RMS.I(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.I);
                   UKF_TP.IE(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.IE);
                   UKF_RMS.IE(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.IE);
                   UKF_TP.P(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.P);
                   UKF_RMS.P(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.P);
                   UKF_TP.PE(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempTP.PE);
                   UKF_RMS.PE(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = mean(tempRMS.PE);                   
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).TP = tempTP;
                    tempVal(clutterLoop,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart).RMS = tempRMS;
                end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            end
            save('loopParam','clutterVersion','totalSim');
            if nSim == totalSim
                save('performance','UKF_TP','UKF_RMS','tempVal');
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart MaxAcc;
            else
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart tempTP tempRMS MaxAcc;
            end 
        end
        disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
    end
    disp(['cluttterLoop:',num2str(clutterLoop)]);    
    clearvars -except MaxAcc;
end
