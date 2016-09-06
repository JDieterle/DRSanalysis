function ok=write_file(filename, data, label, thickness)
labelstr_=label;
if nargin==4
for index=1:numel(thickness)
    labelstr_=[labelstr_, '\t', num2str(thickness(index), '%3.1f'), 'nm'];
end
end

    file=fopen(filename, 'w');
    fprintf(file, labelstr_);
    fprintf(file, '\n');
    fclose(file);
    dlmwrite(filename,data, 'delimiter', '\t', 'newline', 'pc','precision', 6 , '-append');
end