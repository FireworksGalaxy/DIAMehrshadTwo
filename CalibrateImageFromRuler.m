function calibration = CalibrateImageFromRuler(imageFile, rulerLengthMm)
% CalibrateImageFromRuler - Measure image scale from a ruler segment.
%
% Usage:
%   calibration = CalibrateImageFromRuler();
%   calibration = CalibrateImageFromRuler('Image/IMG_1219.jpeg');
%   calibration = CalibrateImageFromRuler('Image/IMG_1219.jpeg', 10);
%
% First zoom/pan to the ruler area, then draw a line over a known ruler
% length in the image. By default the known length is 10 mm. The returned
% pixelsPerMm value can be used as config.calibration_divisor for diameter
% values in mm.

    if nargin < 1 || isempty(imageFile)
        [fileName, folderName] = uigetfile( ...
            {'*.jpg;*.jpeg;*.png;*.tif;*.tiff;*.bmp', 'Image files'; '*.*', 'All files'}, ...
            'Select calibration image');
        if isequal(fileName, 0)
            error('CalibrateImageFromRuler:NoImageSelected', 'No calibration image was selected.');
        end
        imageFile = fullfile(folderName, fileName);
    end

    if nargin < 2 || isempty(rulerLengthMm)
        rulerLengthMm = 10;
    end

    if rulerLengthMm <= 0
        error('CalibrateImageFromRuler:InvalidLength', 'rulerLengthMm must be greater than zero.');
    end

    imageData = imread(imageFile);

    figureHandle = figure('Name', 'Calibration - ruler measurement');
    axesHandle = axes('Parent', figureHandle);
    imshow(imageData, 'Parent', axesHandle);

    title(axesHandle, 'Zoom/pan to the ruler, then press Enter in the Command Window.');
    zoom(figureHandle, 'on');
    pan(figureHandle, 'on');

    fprintf('Zoom/pan to the ruler in the figure window.\n');
    fprintf('When the ruler is clearly visible, press Enter here in the Command Window.\n');
    pause;

    zoom(figureHandle, 'off');
    pan(figureHandle, 'off');
    title(axesHandle, sprintf('Draw a line along %.4g mm on the ruler, then double-click it.', rulerLengthMm));

    rulerLine = drawline(axesHandle, 'Color', 'yellow');
    wait(rulerLine);

    linePosition = rulerLine.Position;
    positionDifference = diff(linePosition, 1, 1);
    distanceInPx = hypot(positionDifference(1), positionDifference(2));

    if distanceInPx == 0
        error('CalibrateImageFromRuler:ZeroLengthLine', 'The calibration line length is zero pixels. Draw a longer line.');
    end

    calibration = struct();
    calibration.imageFile = imageFile;
    calibration.rulerLengthMm = rulerLengthMm;
    calibration.distanceInPx = distanceInPx;
    calibration.pixelsPerMm = distanceInPx / rulerLengthMm;
    calibration.mmPerPixel = rulerLengthMm / distanceInPx;
    calibration.linePosition = linePosition;

    fprintf('Measured %.2f pixels over %.4g mm.\n', distanceInPx, rulerLengthMm);
    fprintf('Pixels per mm: %.6f\n', calibration.pixelsPerMm);
    fprintf('mm per pixel: %.6f\n', calibration.mmPerPixel);
end
