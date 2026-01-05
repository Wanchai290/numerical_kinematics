  function iTk=th_rotz(theta_rad)
    c=cos(theta_rad); s=sin(theta_rad);
    iTk=[[ c,-s, 0,0]
       [ s, c, 0,0]
       [ 0, 0, 1,0]       
       [ 0, 0, 0,1]       
      ]; 
  end
