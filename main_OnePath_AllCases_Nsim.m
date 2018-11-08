close all;
clear;
load('MaxAcc.mat');
clutterVersion = [100000];
for clutterLoop = 1:length(clutterVersion)
    yawDevRateLoopMin = 4.5;
    yawDevRateLoopJump = .5;
    yawDevRateLoopMax = 7;
    maxPings = 150;
    nYaw = (yawDevRateLoopMax/yawDevRateLoopJump)+1;
    nYawStart = yawDevRateLoopMin/yawDevRateLoopJump;
    if clutterLoop == 1
        delete performance.mat
        delete loopParam.mat
        delete DataStore.mat
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
       
       UKF_Ping.C = zeros(maxPings,nYaw-nYawStart);
       UKF_Ping.F = zeros(maxPings,nYaw-nYawStart);
       UKF_Ping.I = zeros(maxPings,nYaw-nYawStart);
       UKF_Ping.P = zeros(maxPings,nYaw-nYawStart);
       UKF_Ping.IE = zeros(maxPings,nYaw-nYawStart);
       UKF_Ping.PE = zeros(maxPings,nYaw-nYawStart);
       save('performance','UKF_TP','UKF_RMS','UKF_Ping');        
    end
    for yawDevRateLoop = yawDevRateLoopMin:yawDevRateLoopJump:yawDevRateLoopMax
        totalSim = 3;        
        if yawDevRateLoop == yawDevRateLoopMin
            DataStore = Init_Data_Store(totalSim,nYaw-nYawStart);
        end
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
        
        tempPing.C = zeros(maxPings,1);
        tempPing.F = zeros(maxPings,1);
        tempPing.I = zeros(maxPings,1);
        tempPing.IE = zeros(maxPings,1);
        tempPing.P = zeros(maxPings,1);
        tempPing.PE = zeros(maxPings,1);
        for nSim = 1:totalSim
            if nSim == 1 && yawDevRateLoop == yawDevRateLoopMin && clutterLoop == 1
                save('loopParam','clutterVersion','totalSim')                   
            else                
                load('loopParam');                
            end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Struct_Param_Path_Init(turnState,turnRate,yawDev,targetStartPoint,tLimit,vicinityLimit)
%             ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],400,100);
            ParPath = Struct_Param_Path_Init(true,yawDevRateLoop,90,[6000,4000],2*maxPings,100);
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
                
                dummyTrackPings = zeros(maxPings,1);
                dummyTrackPings(trackCommInd) = 1;
                
                switch caseLoop
                    case "C" 
                        tempTP.C(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.C(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.C = tempPing.C + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).CW.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).CW.maxPing = impLength;
                    case "F"
                        tempTP.F(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.F(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.F = tempPing.F + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).FM.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).FM.maxPing = impLength;
                    case "I"
                        tempTP.I(nSim) = (trackPoints/impLength)*100;        
                        tempRMS.I(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.I = tempPing.I + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).I.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).I.maxPing = impLength;
                    case "P"
                        tempTP.P(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.P(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.P = tempPing.P + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).P.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).P.maxPing = impLength;
                    case "IE"                                            
                        tempTP.IE(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.IE(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.IE = tempPing.IE + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).IE.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).IE.maxPing = impLength;
                    case "PE"
                        tempTP.PE(nSim) = (trackPoints/(impLength))*100;        
                        tempRMS.PE(nSim) = sqrt(sum(diffBtwTrueExp(trackCommInd).^2)...
                                                /length(trackCommInd));
                                            
                        tempPing.PE = tempPing.PE + dummyTrackPings;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).PE.Pings = trackCommInd;
                        DataStore((yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart,nSim).PE.maxPing = impLength;
                end
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

                UKF_Ping.C(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.C ./ totalSim)*100;
                UKF_Ping.F(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.F ./ totalSim)*100;
                UKF_Ping.I(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.I ./ totalSim)*100;
                UKF_Ping.P(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.P ./ totalSim)*100;
                UKF_Ping.IE(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.IE ./ totalSim)*100;
                UKF_Ping.PE(:,(yawDevRateLoop/yawDevRateLoopJump)+1 - nYawStart) = (tempPing.PE ./ totalSim)*100;
            end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          
            save('loopParam','clutterVersion','totalSim');
            if nSim == totalSim
                save('performance','UKF_TP','UKF_RMS','UKF_Ping','tempVal');
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart MaxAcc maxPings DataStore;
            else
                clearvars -except clutterLoop yawDevRateLoop yawDevRateLoopMin ...
                yawDevRateLoopJump nYawStart tempTP tempRMS tempPing MaxAcc maxPings DataStore;
            end 
        end
        disp(['yawDevRateLoop:', num2str(yawDevRateLoop)]);
    end
    save('DataStore','DataStore')
    disp(['cluttterLoop:',num2str(clutterLoop)]);    
    clearvars -except MaxAcc;
end

function DataStore = Init_Data_Store(totalSim,length)
%     for i = 1:length
        DataStore(length,totalSim).CW.Pings = 0;
        DataStore(length,totalSim).CW.maxPing = 0;
        
        DataStore(length,totalSim).FM.Pings = 0;
        DataStore(length,totalSim).FM.maxPing = 0;
        
        DataStore(length,totalSim).I.Pings = 0;
        DataStore(length,totalSim).I.maxPing = 0;
        
        DataStore(length,totalSim).P.Pings = 0;
        DataStore(length,totalSim).P.maxPing = 0;
        
        DataStore(length,totalSim).IE.Pings = 0;
        DataStore(length,totalSim).IE.maxPing = 0;
        
        DataStore(length,totalSim).PE.Pings = 0;
        DataStore(length,totalSim).PE.maxPing = 0;

end