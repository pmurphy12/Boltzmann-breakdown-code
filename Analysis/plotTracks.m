function [tracks,I] = plotTracks(dataCat,params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    N = 50; %Number of trajectories to plot
    SAVE = true;

    I = randi(length(dataCat{1}.o),[N,1]);
    tracks.x = [];
    tracks.y = [];
    tracks.t = [];
    tracks.id = [];
    
    
    for i = 1:length(dataCat)
        tracks.x = [tracks.x ; dataCat{i}.pos(I,1)];
        tracks.y = [tracks.y ; dataCat{i}.pos(I,2)];
        tracks.t = [tracks.t ; dataCat{i}.t*ones(length(I),1)];
        tracks.id = [tracks.id ; I];
    end
    
    h = figure; hold on; axis([0 params.lengthGx 0 params.lengthGy])
    sz = linspace(1,100,min([params.Nt+1,101]));
    for i = I'
        scatter(tracks.x(tracks.id == i),tracks.y(tracks.id == i),sz)
    end
    
    if SAVE
        save('tracks','tracks','I')
        savefig(['Trajectories of ' num2str(N) ' cells'])
        saveas(h,['Trajectories of ' num2str(N) ' cells.png'],'png')
    end
    
    
end

