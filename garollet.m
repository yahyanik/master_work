function newgen = garollet(ex,Ptedad,options)

ex = ex(:,1);
rollet_wheel = cumsum(ex) / Ptedad;
newgen = zeros(1,Ptedad);
w = 1/Ptedad;
newans = rand * w;
i = 1; 

for ii = 1:Ptedad 
    for j = i:length(rollet_wheel)
        if(newans < rollet_wheel(j)) 
            newgen(ii) = j;
            i = j;
            break;
        end
    end
    newans = newans + w; 
end
   