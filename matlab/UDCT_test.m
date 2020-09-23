clear,clc
I = zeros(8,8);
I(8,8) = 1;
imagePath = fullfile('..', 'image_dir', 'image.jpg');
I_STRUCT = jpeg_read(imagePath);

fun = @(x)x.data .*I_STRUCT.quant_tables{1};
I_DCT = blockproc(I,[8 8],fun);
fun=@(x)idct2(x.data);
I = blockproc(I_DCT,[8 8],fun);

I_spatial = zeros(15,15);
I_spatial(5:12,5:12) = I;

show_cost_dis_color((I_spatial))


%% Prepare for computing DCT bases
k=0:7;
l=0:7;
[k,l]=meshgrid(k,l);

A=0.5*cos(((2.*k+1).*l*pi)/16);
A(1,:)=A(1,:)./sqrt(2);
A=A';

for mode_r = 8
    for mode_c = 8
        DCTbase = A(:,mode_r)*A(:,mode_c)';
        R = conv2(I_spatial, DCTbase, 'symmetric');
    end
end

show_cost_dis_color((R))