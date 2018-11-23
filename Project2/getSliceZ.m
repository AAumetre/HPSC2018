% Get values at slice 'z' for the time step 'timestep'
% /!\ File must be located in directory 'results' and must be named
% 'c_X.dat' where 'X' is the timestep number.
%
function v = getSliceZ(timeStep, z)
    fid = fopen(['c_' num2str(timeStep) '.dat'], 'r');
    N = fread(fid,1,'int32');
    data = fread(fid,N*N*N,'double');
    v = zeros(N,N);
    for a = 1:1:N
        for b = 1:1:N
            v(a,b) = data(b + N*(a-1)+ N*N*z);
        end
    end
    fclose(fid);
end