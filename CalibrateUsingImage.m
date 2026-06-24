function DistanceInPx=CalibrateUsingImage(I)
figure
imshow(I)
d = drawline;
pos = d.Position;
diffPos = diff(pos);
DistanceInPx = hypot(diffPos(1),diffPos(2));
end