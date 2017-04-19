function [DDb] = DDb_eval(x,parameters)

%{
DDb(1,:,:) = [ 0, 0, 0;
               0, 0, 1;
               0, 1, 0];
   
DDb(2,:,:) = [ 0, 0,-1;
               0, 0, 0;
               1, 0, 0];
   
DDb(3,:,:) = [ 0,-1, 0;
               1, 0, 0;
               0, 0, 0];
%}

DDb(1,:,:) = [ 0, 0, 0 ;
               0, 0 -1 ;
               0,-1, 0];
   
DDb(2,:,:) = [ 0, 0, 1 ;
               0, 0, 0 ;
               1, 0, 0];
   
DDb(3,:,:) = [ 0, 1, 0 ;
               1, 0, 0 ;
               0, 0, 0];

end