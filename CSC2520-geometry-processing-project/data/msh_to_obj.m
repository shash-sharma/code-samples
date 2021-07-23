%% Read msh file

filename = 'test_case_3';
full_name = [filename, '.msh'];
data = dlmread(full_name);

%% Process

% test_case_0
Nv = 120;

% test_case_1
Nv = 408;

% test_case_2
Nv = 6000;

% test_case_3
Nv = 738;

Nf = size(data,1) - Nv;

vdata = data(1:Nv,2:4);
fdata = data(Nv+1:end,6:8);

%% Write to obj file

fileID = fopen([filename, '.obj'],'w');
fprintf(fileID,'v %f %f %f\n', vdata');
fprintf(fileID,'f %d %d %d\n', fdata');
fclose(fileID);


