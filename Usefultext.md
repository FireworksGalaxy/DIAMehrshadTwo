#Markdown
fig = findobj('Type', 'figure', 'Name', 'PostProcessFunction - object seperation');
exportgraphics(fig, 'object_separation_HD.png', 'Resolution', 600);
