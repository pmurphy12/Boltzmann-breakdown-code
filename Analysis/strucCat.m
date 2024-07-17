function [dataCat] = strucCat(dataArray)
% concatenate data from simulation repeats

    sz = size(dataArray);

    if sz(1)*sz(2) == max(sz)
        dataCat = dataArray;
        return
    end
    
    dataCat = dataArray(1,:);
    
    for i = 1:sz(2)
        for ii = 2:sz(1)
            dataCat{i}.pos = cat(1,dataCat{i}.pos,dataArray{ii,i}.pos);
            dataCat{i}.o = cat(1,dataCat{i}.o,dataArray{ii,i}.o);
        end
    end
    
end

