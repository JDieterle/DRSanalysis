function ok=write_E_e1_e2_in_VASE_readable_file(data, filename)
title='film';
try
dlmwrite(filename,title,'delimiter','','newline', 'pc');
dlmwrite(filename, 'eV','delimiter','','newline', 'pc', '-append');
dlmwrite(filename, 'e1e2','delimiter', '','newline', 'pc', '-append');
dlmwrite(filename,data,'delimiter', '\t', 'newline', 'pc','precision', 4 , '-append');
catch
    disp(['write_E_e1_e2_in_VASE_readable_file: unable to write file: ', filename])
end

try
 fileserver=evalin('base', 'fileserver');
 fileserver_host=evalin('base', 'fileserver_host');
 eval(['! scp ', filename,' ', fileserver_host, ':', fileserver]);
end
end