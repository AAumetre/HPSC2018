% Build a movie of solution at slice 'z'
%
function v = movieSliceZ(range, z, fixed)
    mov = VideoWriter('movie.avi');
    open(mov);
    minV = 0;
    maxV = 0;
    for i = range
        v = getSliceZ(i, z);
        if i == range(1)
            maxV = max(max(v));
        end
        figure1 = figure('visible','off');
        if fixed == 1
            axes1 = axes('Parent',figure1,'Layer','top','CLim',[minV maxV]);
        else
            axes1 = axes('Parent',figure1,'Layer','top');
        end
        box(axes1,'on');
        hold(axes1,'on');
        image(v,'Parent',axes1,'CDataMapping','scaled');
        colorbar('peer',axes1);
        frame = getframe(figure1);
        writeVideo(mov,frame);
    end
    close(mov);
end