function animate(tmax, plotfunction, filename)
    % given a plotting function (which should take the format
    % plotfunc(i), where i is the frame number, generates
    % a movie file and saves it to filename
    
    fig = figure;
    im = cell(tmax);

    try
        for t = 1:tmax
            plotfunction(t)
            drawnow
            frame = getframe(fig);
            im{t} = frame2im(frame);
        end
    catch exception  %#ok<NASGU> % tab closed
        rethrow(exception);
        disp('Cancelled!')
        return
    end
    close;
    
    % save as gif
    for i = 1:tmax
        [A,map] = rgb2ind(im{i},256);
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
        end
    end
end
