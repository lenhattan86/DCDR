close all;   

legend_str = []; 
     for i=1:10, 
          x = 1:i:(10*i); 
          plot(x); hold on; 
          legend_str = [legend_str; {num2str(i)}];
     end
     columnlegend(3, legend_str, 'Location', 'NorthWest');