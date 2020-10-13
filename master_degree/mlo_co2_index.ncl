;------------------------------------
  begin
  diri="/mnt/f/share/CO2/"
  params=systemfunc("ls "+diri+"MLO_CO2_annual.txt")
  ncols =3
  nrows=58
  co2_r=asciiread(params, (/nrows,ncols/), "float")
  co2=co2_r(0:,1)
  co2_std=dim_standardize_n(co2, 1, 0)
  co2_std!0="time"
  co2_std&time=ispan(1959,2016,1)
  res=True

  wks=gsn_open_wks("png", "/mnt/f/share/era/data/co2_index")
  plot=gsn_csm_xy(wks,co2_std&time , co2_std, res)
 
  opt="/mnt/f/share/era/data/mlo_co2_annual_index_1959_2016.ascii"    ;;;输出强度指数
  asciiwrite(opt, sprintf ("%6.2f",co2_std))  ;;;输出强度指数
  
  
  end 