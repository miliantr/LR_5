function out = read_data(file_name, cols)
    f_id = fopen(file_name,'r');
    formatSpes = "%f";
    out = fscanf(f_id, formatSpes, [cols, inf]);
    fclose(f_id);
    out = out';
end