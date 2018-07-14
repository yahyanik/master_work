function stop = outfun(x,optimValues,state) 

stop = false;
StopFlag = getappdata(0,'FMINCONsStopFlag');
if StopFlag 
stop = true; 
end
if isnan(x) 
    stop = true; 
end
end