function newans  = crossoverfunc(parents,options,GenomeLength,~,~,thisPopulation)
nKids = length(parents)/2;
newans = zeros(nKids,GenomeLength);
for i=1:2:nKids
    r1 = parents(i);
    j = i + 1;
    r2 = parents(j);
    for j = 1:GenomeLength
        if(rand > 0.5)
            newans(i,j) = thisPopulation(r1,j);
        else
            newans(i,j) = thisPopulation(r2,j);
        end
    end
end
