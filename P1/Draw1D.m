function Draw1D( ER, E, H, dz )

   za=[0:length(E)-1]*dz;
   basecolor = [1 1 1];
   Color = 1./ER;

   fill(0,0,'-w');   hold on;
   i = 1;
   count = 0;
   prev = 0;
   while i < length(ER)
    i = i + 1;
    
    if(prev == 1)
      prev = ER(i);
      continue;
    end
    
    if(prev == ER(i))
      count = count + 1;
      
    else
      xstart = (i-count)*dz;
      xend = xstart + count*dz;
      
      x = [ xstart xend xend xstart xstart ];
      y = [ -1.5 -1.5 1.5 1.5 -1.5 ];
      fill(x,y,[Color(i-1) Color(i-1) Color(i-1)]);
      
      count = 1;
      prev = ER(i);
    end
   end

   plot(za, E, '-b'); 
   plot(za, H, '-r'); 
   
   hold off;
end

