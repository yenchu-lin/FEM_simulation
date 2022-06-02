% Write image to FPI (binary) file
% For EyeRIS experiment
% filename is the directory/filename.fpi that you want to save to
% img is the matrix [m x n] containing the image you would like to save
% fid2 = fopen(filename, 'w');
% fwrite(fid2, img, 'float');
% fclose(fid2);
% (K is the img here)

letter = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

for ii = 1:length(letter)
    
    load_filename = [letter{ii},'.mat'];
    S = load(load_filename);
    fields = fieldnames(S);
    
    save_filename = ['letter',letter{ii}];
    fid = fopen(save_filename, 'w');
    fwrite(fid, double(S.(fields{1})), 'float');
    fclose(fid);
    
end

%% create blank
load_filename = 'blank.mat';
S = load(load_filename);
fields = fieldnames(S);

save_filename = 'blank';
fid = fopen(save_filename, 'w');
fwrite(fid, double(S.(fields{1})), 'float');
fclose(fid);
