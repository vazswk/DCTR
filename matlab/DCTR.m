function F = DCTR(I_STRUCT, QF)

%% Set parameters
% number of histogram bins
T = 4;

% compute quantization step based on quality factor
if QF<50,
    q = min(8 * (50 / QF), 100);
else
    q = max(8 * (2 - (QF/50)), 0.2);
end

%% Prepare for computing DCT bases
k=0:7;
l=0:7;
[k,l]=meshgrid(k,l);

A=0.5*cos(((2.*k+1).*l*pi)/16);
A(1,:)=A(1,:)./sqrt(2);
A=A';

%% Compute DCTR locations to be merged
% coordinates<9 : <9 的位置输出1，否则输出0；
% 判断x的每一行是否不包含0元素，使用all(x,2)：输出1代表该行不包含0，输出0代表该行包含0

mergedCoordinates = cell(25, 1);
for i=1:5
    for j=1:5
        coordinates = [i,j; i,10-j; 10-i,j; 10-i,10-j];
        coordinates = coordinates(all(coordinates<9, 2), :);  %  结合上一句，实现了，输出coordinates中数值小于9的所有行
        mergedCoordinates{(i-1)*5 + j} = unique(coordinates, 'rows'); % 如果coordinates中有重复的行，那么只输出不同行的结果
    end
end

%% Decompress to spatial domain
fun = @(x)x.data .*I_STRUCT.quant_tables{1};
I_spatial = blockproc(I_STRUCT.coef_arrays{1},[8 8],fun);
fun=@(x)idct2(x.data);
I_spatial = blockproc(I_spatial,[8 8],fun)+128;

%% Compute features
modeFeaDim = numel(mergedCoordinates)*(T+1);% 125
F = zeros(1, 64*modeFeaDim, 'single');
for mode_r = 1:8
    for mode_c = 1:8
        modeIndex = (mode_r-1)*8 + mode_c;
        % Get DCT base for current mode
        DCTbase = A(:,mode_r)*A(:,mode_c)';
        
        % Obtain DCT residual R by convolution between image in spatial domain and the current DCT base
        R = conv2(I_spatial-128, DCTbase, 'valid');
                
        % Quantization, rounding, absolute value, thresholding
        R = abs(round(R / q));       % 量化
        R(R > T) = T;    % 截断
        
        % Core of the feature extraction
        for merged_index=1:numel(mergedCoordinates)
            f_merged = zeros(1, T+1, 'single');
            for coord_index = 1:size(mergedCoordinates{merged_index}, 1);
                r_shift = mergedCoordinates{merged_index}(coord_index, 1);
                c_shift = mergedCoordinates{merged_index}(coord_index, 2);
                R_sub = R(r_shift:8:end, c_shift:8:end);
                f_merged = f_merged + hist(R_sub(:), 0:T);
            end
            F_index_from = (modeIndex-1)*modeFeaDim+(merged_index-1)*(T+1)+1;
            F_index_to = (modeIndex-1)*modeFeaDim+(merged_index-1)*(T+1)+T+1;
            F(F_index_from:F_index_to) = f_merged / sum(f_merged); % 归一化
        end
    end
end

end

