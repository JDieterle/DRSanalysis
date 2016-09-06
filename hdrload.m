function [data, data2] = hdrload(file)
try
data=[];
data2=load(file);
catch
    disp('usage of standard hdrload.m file can solve some problems')
end
end

