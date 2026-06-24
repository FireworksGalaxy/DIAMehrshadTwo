function new_image=MaskOverlappedBubblesBeforeWatershed(I3,masked,config)
% Step 1: Extract the background color (use mean color of unmasked areas)

background_mask = ~masked;  % Invert mask to get background region

red_channel = I3(:,:,1);
green_channel = I3(:,:,2);
blue_channel = I3(:,:,3);

bg_r = mean(red_channel(background_mask)); 
bg_g = mean(green_channel(background_mask));
bg_b = mean(blue_channel(background_mask));


% Step 2: Create a new image with the same background color
[rows, cols, ~] = size(I3);
new_background = uint8(cat(3, ...
    ones(rows, cols) * bg_r, ...
    ones(rows, cols) * bg_g, ...
    ones(rows, cols) * bg_b));

% Step 3: Extract the two bubbles
bubbles = I3;
bubbles(~masked) = 0;  % Keep only the bubbles, set background to black

% Step 4: Place the bubbles onto the new background
new_image = new_background;
new_image(masked) = bubbles(masked);  % Overlay bubbles onto background

% Display the result
if config.show_figures
    figure('Name','First Image After Masking and Before Watershedding');
    imshow(new_image);
end
end
