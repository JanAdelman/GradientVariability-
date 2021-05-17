I = imread('plot1.png');
close all
imshow(I, 'XData', [0 120], 'YData', [0.16 0])
axis on
axis normal
set(gca, 'YDir', 'normal')
datacursormode
