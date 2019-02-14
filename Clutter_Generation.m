 function ClutterStruct = Clutter_Generation(TarMaxRho,clutterDensity,impInstVec)
    
%     areaSquare = (2*TarMaxRho)^2;
%     ptsGen = ceil(areaSquare*clutterDensity);
%     n = size(impInstVec,1);
%     xClutter = -TarMaxRho + 2*TarMaxRho*rand(ptsGen,n);          %This is uniform distributiom
%     yClutter = -TarMaxRho + 2*TarMaxRho*rand(ptsGen,n);          %This is uniform distributiom
%     
%     for i = n:-1:1
%         [thetaClutter,rhoClutter] = cart2pol(xClutter(:,i),yClutter(:,i));
%         indx = find(rhoClutter<=TarMaxRho);
%         ClutterStruct(i).thetaClutter = thetaClutter(indx,:);
%         ClutterStruct(i).rhoClutter = rhoClutter(indx,:);
%         noClutterPts = length(indx);        
%         ClutterStruct(i).rangeRateClutter = 5*randn(noClutterPts,1);
%     end
    area = pi*TarMaxRho^2;
    points = ceil(area*clutterDensity);
    n = size(impInstVec,1);
    thetaClutter = 2*pi*rand(points,n);                                %This is uniform distributiom
    rhoClutter = TarMaxRho*rand(points,n);                               %This is uniform distributiom
    rangeRateClutter = 5*randn(points,n);
    for i = n:-1:1       
            ClutterStruct(i).thetaClutter = thetaClutter(:,i);
            ClutterStruct(i).rhoClutter = rhoClutter(:,i);              
            ClutterStruct(i).rangeRateClutter = rangeRateClutter(:,i);
    end