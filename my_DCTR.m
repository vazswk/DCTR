function my_DCTR(image_path, feature_path, QF)

QF = uint32(QF);
files = dir([image_path '\*.jpg']);
file_num = length(files);
F = single( zeros(file_num,8000));
names = cell(file_num,1);

for w =1:length(files)
    tic
    jpegfilename = [image_path '\' files(w).name]; % jpegfilename是string数据类型
    ImageSet = {jpegfilename}; % ImageSet是cell数据类型
    f = DCTR(ImageSet, QF); % f是single数据类型
    F(w,:) = f(:); % F是single数据类型
    % names{w} = [num2str(w) '.jpg'];
    names{w} = files(w).name;
    toc
end
save(feature_path,'F','names');
% save(feature_path,'F','names','-v7.3');
