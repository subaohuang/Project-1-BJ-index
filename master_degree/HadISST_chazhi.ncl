begin
;;;读入
	diri="/mnt/f/share/ersst/"
	diri_pi="/mnt/f/share/era/"
	f0=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
	date=f0->time
	lon=f0->longitude
	lat=f0->latitude

	time=cd_calendar(date, 1)
	timeind=ind(time.le.201212.and.time.ge.197901)
	sst_r=lonFlip(f0->sst(timeind,:,:))
;;;
 
  
;;;插值
	xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

	 sst_int=linint2_Wrap(sst_r&longitude, sst_r&latitude(::-1), sst_r(:,::-1,:), False, xo, yo, 0)

	 printVarSummary(sst_int)
    sst=sst_int

    sst!0="time"
    sst!1="lat"
    sst!2="lon"
    sst&lat@units="degree_north"
    sst&lon@units="degree_east"
     printVarSummary(sst)
;;;
     

;****************************************************************************
;随时间变化的量，nc文件
    diro=diri
    filo1="HadISST_sst_chazhi_2012.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "HadISST_sst插值的结果。插值到89-89,1到359"    
    fAtt@creation_date = "2017.10.25"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=getVarDimNames(sst)
    dimsize=dimsizes(sst)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "lat",typeof(sst&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(sst&lon),(/"lon"/))
    filevardef(fout, "time", typeof(sst&time),(/"time"/))
    filevardef(fout, "sst", typeof(sst),(/"time","lat","lon"/))    


    filevarattdef(fout, "sst", sst)
    filevarattdef(fout, "time", sst&time)
    filevarattdef(fout, "lat", sst&lat)
    filevarattdef(fout, "lon", sst&lon)

    setfileoption(fout,"DefineMode",False)

    fout->lat=(/sst&lat/)
    fout->lon=(/sst&lon/)
    fout->time=(/sst&time/)
    fout->sst=(/sst/)
;**********************************************************************   
end
