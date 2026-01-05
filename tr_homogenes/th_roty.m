 function iTk=th_roty(theta_rad)
    c=cos(theta_rad); s=sin(theta_rad);
    iTk=[[ c, 0, s,0]
       [ 0, 1, 0,0]
       [-s, 0, c,0]       
       [ 0, 0, 0,1]       
      ]; 
  end
