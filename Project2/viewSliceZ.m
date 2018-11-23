% Same as getSliceZ but draws the solution
%
function viewSliceZ(timeStep, z)
    fid = fopen(['results/c_' num2str(timeStep) '.dat'], 'r');
    N = fread(fid,1,'int32');
    data = fread(fid,N*N*N,'double');
    v = zeros(N,N);
    for a = 1:1:N
        for b = 1:1:N
            v(a,b) = data(b + N*(a-1)+ N*N*z);
        end
    end
    fclose(fid);
    figure1 = figure;
    axes1 = axes('Parent',figure1,'Layer','top');
    box(axes1,'on');
    hold(axes1,'on');
    image(v,'Parent',axes1,'CDataMapping','scaled');
    colorbar('peer',axes1);
end