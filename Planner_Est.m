function transWaveEst = Planner_Est(FilterPredict,ParPath,caseLoop,MVec,inst)
  
    switch caseLoop
        case "C"
           transWaveEst = caseLoop;
        case "F"
            transWaveEst = caseLoop;
        case "I"
            estTargetDirTheta = Est_Target_Direction(FilterPredict.XkMinus,ParPath.selfX,ParPath.selfY);
            srcWaves = ['F','C'];
            chanceNodeVec = zeros(1,size(srcWaves,2));
            for i = 1: length(srcWaves)
                chanceNodeVec(i) = Est_Eval(estTargetDirTheta,srcWaves(i));
            end
            transWaveEst = string(srcWaves(chanceNodeVec == max(chanceNodeVec))); 
        case "IE"
            if inst < 2
                transWaveEst = Planner_Est(FilterPredict,ParPath,"I",MVec,inst);
            else
                if (MVec(inst-1) > 7)% && (MVec(inst-2,MVec(inst-2,4)) > 7)                
                    transWaveEst = "C";                   
                else
                    transWaveEst = Planner_Est(FilterPredict,ParPath,"I",MVec,inst);
                end
            end
        case "PE"
            if inst ==1
                transWaveEst = Planner_Est(FilterPredict,ParPath,"I",MVec,inst);
            elseif inst ==2
                if (MVec(inst-1) > 7)            
                    transWaveEst = "C";                   
                else
                    transWaveEst = Planner_Est(FilterPredict,ParPath,"I",MVec,inst);
                end                
            else
                if (MVec(inst-2) > 7)      
                    transWaveEst = "C";                   
                else
                    transWaveEst = Planner_Est(FilterPredict,ParPath,"I",MVec,inst);
                end
            end
    end     
end


function estTargetDirTheta = Est_Target_Direction(XkMinus,selfX,selfY)

    velDirSlope         = XkMinus(4)/XkMinus(3);
    slopeWithOrigin     = (XkMinus(2)- selfY)./(XkMinus(1)- selfX);
    estTargetDirTheta   = atan(abs((velDirSlope-slopeWithOrigin)...
                         ./(1+ (velDirSlope.*slopeWithOrigin))));
end

function score = Est_Eval(estTargetDirTheta,srcWave)

    switch srcWave
        case 'F'
            fm = 1;
            cw = 0;
        case 'C'
            fm = 0;
            cw = 1;
    end
    score = ((cos(estTargetDirTheta))^2) * cw + ((sin(estTargetDirTheta))^2)* fm;    
end
