begin
 diri1="/mnt/f/share/era/data/"
 fils1=systemfunc("ls "+diri1+"era_u_sst_slp.nc")
 f1=addfile(fils1, "r")
 time=f1->time
 lon=f1->longitude
 lat=f1->latitude
 date=cd_calendar(time, 1)

 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201612)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 lonind1=ind(lon.eq.190)
 lonind2=ind(lon.eq.240)

 sst_r=f1->sst(tind1:tind2,latind1:latind2,lonind1:lonind2)
 
 sst=sst_r*sst_r@scale_factor+sst_r@add_offset-273.15
 copy_VarMeta(sst_r, sst)
 ssta_de=dtrend_n(rmMonAnnCycLLT(sst(latitude|:,longitude|:,time|:)),True,2)
 nino3_4=dim_avg_n_Wrap(ssta_de, (/0,1/))
 print(nino3_4)
 opt=diri1+"era_nino3_4.ascii"
 asciiwrite(opt, sprintf ("%5.2f,", nino3_4))
end