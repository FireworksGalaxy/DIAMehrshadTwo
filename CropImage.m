addpath 'C:\Users\20233204\OneDrive - TU Eindhoven\Desktop\phd\Report\Report 23\DIA\Second Start';
img = imread('0.6LowResolution.jpg'); % Read the image
imshow(img); % Display the image
cropped_img = imcrop(img); % Select and crop using the mouse
imshow(cropped_img); % Display the cropped image
ims