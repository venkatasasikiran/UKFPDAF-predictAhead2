function [FilterEstimateUKF,transWave,Measurements] = UKF_PDAF3(ParPath,ParGen,Timings,Measurements,...
                                                                noiseAcc,ClutterStruct,ParStdDev)
    
    n = size(Timings.impInstVec,1);
    extImpInstVec = [1;Timings.impInstVec];
    stateSize = 4;
    transWave = string(zeros(n,1));
    FilterEstimateUKF = Struct_Param_Filter_Init(ParPath,stateSize,n);    
    ParPDAF = Struct_Param_PDAF_init;
    
    Measurements.MeasFinal.Pos = zeros(n,2);
    Measurements.MeasFinal.rangeRate = nan(n,1);
    MVec = zeros(n,1);
    
    for i = 1:n
        if i == 1
            T = (extImpInstVec(i+1) - extImpInstVec(i))/ParGen.f;
            FilterMat = Struct_Mat_Filter(T,ParStdDev,noiseAcc);
            FilterEstimateUKF = UKF_Estimate(FilterEstimateUKF,FilterMat,i);         
            transWave(i,:) = Planner_Est(FilterEstimateUKF,ParPath,"I",MVec,i);
        end        
        [FilterEstimateUKF,z,M] = UKF_Update(FilterEstimateUKF,FilterMat,Measurements, ...
                                   ClutterStruct,transWave(i),ParPDAF,i);
        MVec(i) = M;        
        switch transWave(i)
            case "C"
                Measurements.MeasFinal.Pos(i,:) = z(1:2)';
                Measurements.MeasFinal.rangeRate(i) = z(3);
            case "F"
                Measurements.MeasFinal.Pos(i,:) = z';
        end
        
        if i < n
            T2 = (extImpInstVec(i+2) - extImpInstVec(i+1))/ParGen.f;
            FilterMat = Struct_Mat_Filter(T2,ParStdDev,noiseAcc);
            FilterEstimateUKF = UKF_Estimate(FilterEstimateUKF,FilterMat,i+1);
            if i == 1            
                transWave(i+1,:) = Planner_Est(FilterEstimateUKF,ParPath,"I",MVec,i+1);
            end
        end
        
        if i < n-1
            T3 = (extImpInstVec(i+3) - extImpInstVec(i+2))/ParGen.f;
            FilterMat2 = Struct_Mat_Filter(T3,ParStdDev,noiseAcc);
            FilterEstimateUKF2 = UKF_Estimate2(FilterEstimateUKF,FilterMat2,i+2);        
            transWave(i+2,:) = Planner_Est(FilterEstimateUKF2,ParPath,"I",MVec,i+2);
        end       
        
    end



end
function FilterMat = Struct_Mat_Filter(T,ParStdDev,noiseAcc)

    FilterMat.F = [1 0 T 0;...
                   0 1 0 T;...
                   0 0 1 0;...
                   0 0 0 1];
    FilterMat.L = [0.5*(T)^2 0                 ;...     
                         0 0.5*(T)^2 ;...
                         T 0                 ;...
                         0 T];  
    qX = noiseAcc(1)*1;
    qY = noiseAcc(2)*1;
    Qxy = [qX^2 0;
           0    qY^2];
    FilterMat.Q = FilterMat.L*Qxy*FilterMat.L';
    
    FilterMat.R.Cw = [ParStdDev.Range.cw^2             0               0;  
                               0               ParStdDev.theta^2 0;
                               0                           0               ParStdDev.RangeRate.cw^2];
    FilterMat.R.Fm = [ParStdDev.Range.fm^2             0;  
                               0              ParStdDev.theta^2];
end


function ParPDAF = Struct_Param_PDAF_init()
    ParPDAF.gateGamma = 16;  
    ParPDAF.Pd = 1;
    ParPDAF.Pg.fm = 0.9997;
    ParPDAF.Pg.cw = 0.9989;
    ParPDAF.clutterDensity = 1/100000;
end