
;基于era数据的，定义各强度指数的结果，进行对比。
function gain_data(data1)
 begin
 
 rc1=regCoef(data1&time,data1)

 data_fuc=new((/4,dimsizes(data1&time)/),double)
 data_fuc(0,:) = data1(:)
 data_fuc(1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(2,:) =0.0
 data_fuc(3,:) =rc1
return(data_fuc)
end 
;********************************************************
function clm_ano(x[*]:numeric)  ;去除一维时间序列的年循环
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   xAnom=x
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
   end do
   return (xAnom)
end
;*************************************************************
begin

diri="/mnt/f/share/era/"
;********************************************************
;1：定义的500hpa上msf累加结果，以及采用周天军150E到150W的结果（STRF）
;*******************************************************
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 level=f1->level(:)
 lon=f1->lon(:)
 lonind1=ind(lon.eq.150)
 lonind2=ind(lon.eq.210)
 msf_time=f1->msf_time

 STRF=dim_avg_n_Wrap(msf_time(:,lonind1:lonind2,:),(/0,1/))
 printVarSummary(STRF)
 msf_time_new=where(msf_time.gt.0, msf_time, msf_time@_FillValue)
 copy_VarMeta(msf_time, msf_time_new)
 intesion=dim_sum_n_Wrap(msf_time_new(ind(msf_time_new&level.eq.500),ind(lon.eq.140):ind(lon.eq.260),:), 0)  
 printVarSummary(intesion)

 intesion_an=runave_n_Wrap(dim_rmvmean(month_to_annual(intesion(0:455), 1))/10^10,3,0,0)
 intesion_an!0="time"
 intesion_an&time=ispan(1979,2016,1)
 data=gain_data(intesion)
 data_an=gain_data(intesion_an)

 STRF_an=doubletofloat(runave_n_Wrap(dim_rmvmean(month_to_annual(STRF(0:455), 1))/10^10,3,0,0))
 STRF_an!0="time"
 STRF_an&time=ispan(1979,2016,1)

 STRF_data=gain_data(STRF_an)                                     



 
 

 delete(lon)
 delete(level)
 
;*****************************************************************
  ;2 slp定义的强度变化,取5°S-5°N 160-80W 和达尔文区5S-5N 80-160E的差值。
;**********************************************************************
 fils2=systemfunc("ls "+diri+"era_slp_80E-70w.nc")
 f2=addfile(fils2, "r")

 time=f2->time(:)
 lon=f2->longitude(:)
 lat=f2->latitude(:)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 date=cd_calendar(time, 1)
 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201612)

 msl=f2->msl(tind1:tind2,latind1:latind2,:)
 m=getvaratts(msl)
 print(m)
  if (any(getvaratts(msl).eq."scale_factor")) then 
      msl_new=dble2flt(msl*msl@scale_factor+msl@add_offset)
  end if
  copy_VarCoords(msl, msl_new)
  delete(msl)
 msl_tahiti=dim_avg_n_Wrap(msl_new(latitude|:,longitude|ind(lon.eq.200):ind(lon.eq.280),time|:),(/0,1/))
 msl_daerwen=dim_avg_n_Wrap(msl_new(latitude|:,longitude|ind(lon.eq.80):ind(lon.eq.160),time|:),(/0,1/))
 slp=(msl_tahiti - msl_daerwen)/100

 slp!0="time"
 slp&time=ispan(1,456,1)

 intesion_slp=runave_n_Wrap(dim_rmvmean(month_to_annual(slp, 1)),3,0,0)  ;计算年份，并计算anomoly，再滑动平均
 intesion_slp!0="time"
 intesion_slp&time=ispan(1979,2016,1)
   
 slp_data=gain_data(intesion_slp)
 
delete(time)
delete(lon)
delete(lat)
delete(latind2)
delete(latind1)
delete(date)
delete(m)
;*****************************************************************
  ;3 ud定义的强度，基于平均的结果，150E到210度E，150到300和1000到750的平均结果
;**********************************************************************
 fils3=systemfunc("ls "+diri+"msf_time_era.nc")
 f3=addfile(fils3, "r")
 level=f3->level
 ud_time=f3->ud_time 
 intension_ud=dim_avg_n_Wrap(ud_time(ind(level.eq.300):ind(level.eq.150),lonind1:lonind2,:),(/0,1/))-dim_avg_n_Wrap(ud_time(ind(level.eq.1000):ind(level.eq.750),lonind1:lonind2,:),(/0,1/))

 intension_ud!0="time"
 intension_ud&time=ispan(1,457,1)

 intension_ud_an=runave_n_Wrap(dim_rmvmean(month_to_annual(intension_ud(0:455), 1)),3,0,0)
 intension_ud_an!0="time"
 intension_ud_an&time=ispan(1979,2016,1)
 ud_data=gain_data(intension_ud_an)     
 delete(level)
;*****************************************************************
  ;4 -U10定义的强度，基于平均的结果，150E到210度E，150到300和1000到750的平均结果
;**********************************************************************

 fils4=systemfunc("ls "+diri+"era_u_sst_slp.nc")
 f4=addfile(fils4, "r")

 time=f4->time(:)
 lon=f4->longitude(:)
 lat=f4->latitude(:)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 lonind1=ind(lon.eq.150)
 lonind2=ind(lon.eq.210)
 
 date=cd_calendar(time, 1)
 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201701)
 
 u10=f4->u10(tind1:tind2,latind1:latind2,lonind1:lonind2)

  if (any(getvaratts(u10).eq."scale_factor")) then 
      u10_new=dble2flt(u10*u10@scale_factor+u10@add_offset)
  end if
  copy_VarCoords(u10, u10_new)
  delete(u10)
 u10_mean=dim_avg_n_Wrap(u10_new, (/1,2/))

 intension_u10_an=-1*runave_n_Wrap(dim_rmvmean(month_to_annual(u10_mean(0:455), 1)),3,0,0)
 intension_u10_an!0="time"
 intension_u10_an&time=ispan(1979,2016,1)
 u10_data=gain_data(intension_u10_an)  
 delete(time)
;**************************************************************
;去趋势：
;****************************************************************

    ; intesion_slp_de=dtrend_n(intesion_slp,True,0)
    ; intension_u10_an_de=dtrend_n(intension_u10_an,True,0)
    ; STRF_an_de=dtrend_n(STRF_an,True,0)
    ; intesion_an_de=dtrend_n(intesion_an,True,0)
    ; intension_ud_an_de=dtrend_n(intension_ud_an,True,0)
    ; copy_VarCoords(intesion_slp,intesion_slp_de)
    ; copy_VarCoords(intension_u10_an,intension_u10_an_de)
    ; copy_VarCoords(STRF_an,STRF_an_de)
    ; copy_VarCoords(intesion_an,intesion_an_de)
    ; copy_VarCoords(intension_ud_an,intension_ud_an_de)
    
    ; data_slp_de=gain_data(intesion_slp_de)
    ; data_u10_de=gain_data(intension_u10_an_de)
    ; data_STRF_de=gain_data(STRF_an_de)
    ; data_intesion_de=gain_data(intesion_an_de)
    ; data_ud_de=gain_data(intension_ud_an_de)


;*************************************************************************
 ;画图
;**************************************************************************8
  
;画年际变化：
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString="(a)dSLP"
    ; resxy@gsnRightString=""+slp_data(3,1)+""
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01
    ; resxy@trXMaxF=dimsizes(intesion&time)
    ; resxy@tmXBValues=ispan(13,dimsizes(intesion&time),60)
    ; resxy@tmXBLabels=ispan(1980,2015,5)

    resxy@tmXBMode="Explicit"
    resxy@trXMaxF=2016
    resxy@trXMinF=1979
    resxy@tmXBValues=ispan(1979,2016,5)
    resxy@tmXBLabels=ispan(1979,2016,5)
    resxy@tmXBMinorOn=True
    resxy@tmXBMinorValues=ispan(1979,2016,1)

    resxy@tmYLMinorOn=False

    resxy@vpHeightF=0.13
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.95

    resxy@xyLineThicknesses = (/ 3.0,1.5,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change",  "Linear trend"/)     ;加线的标注

    resxy@xyDashPatterns    =(/0,0,2/)
    
    
    resxy@tmXBLabelFontAspectF=1.0
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    ; resxy@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

    resxy@tmXBMinorLengthF=0.003
    resxy@tmXBMinorOutwardLengthF=0.003

    resxy@tmYLMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    resxy@tmYLLabelFont=30
    resxy@tmYLLabelFontHeightF=0.008
    resxy@tmXTOn=False
    resxy@tmYROn=False
    ; resxy@tiYAxisString = title

    ; resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    ; resxy@pmLegendSide           = "Top"               ; Change location of 
    ; resxy@pmLegendParallelPosF   = 0.85                  ; move units right
    ; resxy@pmLegendOrthogonalPosF = -0.18                ; move units down
    ; resxy@pmLegendWidthF         = 0.10                ; Change width and
    ; resxy@pmLegendHeightF        = 0.1               ; height of legend.
    ; resxy@lgPerimOn              = False               ; turn off box around
    ; resxy@lgLabelFontHeightF     = .008                 ; label font height

    ; resxy@xyExplicitLegendLabels = (/"Anual change", "linear trend"+"("+sprintf("%5.3f",rc1(1)*120)+"~S~o~N~/decade"+")","average"/)


   wks=gsn_open_wks("png", ""+diri+"walker_intension")
   wks1=gsn_open_wks("png", ""+diri+"msf_intension_500hpa_era_an")
   wks2=gsn_open_wks("png", ""+diri+"msf_intension_era_strf_an")
   wks3=gsn_open_wks("png", ""+diri+"intesion_slp_an")
   wks4=gsn_open_wks("png", ""+diri+"intesion_ud_an")
   plot=new(5,graphic)

   plot(0)=gsn_csm_xy(wks,intesion_slp&time,slp_data(0:2,:),resxy)
    resxy@gsnLeftString="(b)-U10"
    ; resxy@gsnRightString=""+u10_data(3,1)+""
      resxy@vpYF=0.77
      
   plot(1)=gsn_csm_xy(wks,intension_u10_an&time,u10_data(0:2,:),resxy)
    resxy@gsnLeftString="(c)STRF"
    ; resxy@gsnRightString=""+STRF_data(3,1)+""
      resxy@vpYF=0.59
      resxy@tmYLPrecision=1
   plot(2)=gsn_csm_xy(wks,STRF_an&time,STRF_data(0:2,:),resxy)
       delete(resxy@tmYLPrecision)
    resxy@gsnLeftString="(d)msf_500"
    ; resxy@gsnRightString=""+data_an(3,1)+""
      resxy@vpYF=0.41
   plot(3)=gsn_csm_xy(wks,intesion_an&time,data_an(0:2,:),resxy)
   resxy@gsnLeftString="(e)ud"
   ; resxy@gsnRightString=""+ud_data(3,1)+""
    resxy@vpYF=0.23
    resxy@tmYLPrecision=1
   plot(4)=gsn_csm_xy(wks,intension_ud_an&time,ud_data(0:2,:),resxy)
      
    
    ; gsn_panel(wks, plot(0:2), (/3,1/), pres)

    enino=(/1982,1983,1986,1987,1991,1992,1993,1994,1997,1998/)
    print(intesion_an&time)
    ilon=new(dimsizes(enino),integer)
    do i=0,dimsizes(enino)-1
    ilon(i)=ind(intesion_an&time.eq.enino(i))
  end do
    itime=intesion_an&time(ilon)
   print(itime)
    ilev=new((/5,dimsizes(enino)/),double)
    ilev(0,:)=intesion_slp(ilon)
    ilev(1,:)=intension_u10_an(ilon)
    ilev(2,:)=STRF_an(ilon)
    printVarSummary(intesion_an)
    ilev(3,:)=intesion_an(ilon)
    ilev(4,:)=intension_ud_an(ilon)
   
   
   polyres                = True
   polyres@gsMarkerIndex     = 16          ; polymarker style
   polyres@gsMarkerSizeF     = 5.         ; polymarker size
   polyres@gsMarkerColor     ="red"
   polyres@tfPolyDrawOrder="PreDraw"
   dum1=new(5,graphic)
   do i=0,4
    dum1(i)=gsn_add_polymarker(wks, plot(i), itime, ilev(i,:), polyres)
   end do 

    pres=True
    maximize_output(wks,pres)  ; calls draw and frame for you 
    

  ;********************************************************************
   ;  wks2=gsn_open_wks("png", ""+diri+"walker_intension_detrend")
   ; plot_de=new(5,graphic)  
   ; plot_de(0)=gsn_csm_xy(wks,intesion_slp_de&time,data_slp_de(0:2,:),resxy)
   ;  resxy@gsnLeftString="(b)-U10"
   ;  ; resxy@gsnRightString=""+u10_data(3,1)+""
   ;    resxy@vpYF=0.77
      
   ; plot_de(1)=gsn_csm_xy(wks2,intension_u10_an&time,data_u10_de(0:2,:),resxy)
   ;  resxy@gsnLeftString="(c)STRF"
   ;  ; resxy@gsnRightString=""+STRF_data(3,1)+""
   ;    resxy@vpYF=0.59
   ;    resxy@tmYLPrecision=1
   ; plot_de(2)=gsn_csm_xy(wks2,STRF_an&time,data_STRF_de(0:2,:),resxy)
   ;     delete(resxy@tmYLPrecision)
   ;  resxy@gsnLeftString="(d)msf_500"
   ;  ; resxy@gsnRightString=""+data_an(3,1)+""
   ;    resxy@vpYF=0.41
   ; plot_de(3)=gsn_csm_xy(wks2,intesion_an&time,data_intesion_de(0:2,:),resxy)
   ; resxy@gsnLeftString="(e)ud"
   ; ; resxy@gsnRightString=""+ud_data(3,1)+""
   ;  resxy@vpYF=0.23
   ;  resxy@tmYLPrecision=1
   ; plot_de(4)=gsn_csm_xy(wks2,intension_ud_an&time,data_ud_de(0:2,:),resxy)
      
    
   ;  ; gsn_panel(wks2, plot_de(0:2), (/3,1/), pres)

   
   ;  ilev_de=new((/5,dimsizes(enino)/),double)
   ;  ilev_de(0,:)=intesion_slp_de(ilon)
   ;  ilev_de(1,:)=intension_u10_an_de(ilon)
   ;  ilev_de(2,:)=STRF_an_de(ilon)
   ;  printVarSummary(intesion_an)
   ;  ilev_de(3,:)=intesion_an_de(ilon)
   ;  ilev_de(4,:)=intension_ud_an_de(ilon)

  
   ; polyres_de_de                   = True
   ; polyres_de_de@gsMarkerIndex     = 16          ; polymarker style
   ; polyres_de_de@gsMarkerSizeF     = 5.         ; polymarker size
   ; polyres_de_de@gsMarkerColor     ="red"
   ; polyres_de_de@tfPolyDrawOrder="PreDraw"
   ; dum2=new(5,graphic)
   ; do i=0,4
   ;  dum2(i)=gsn_add_polymarker(wks2, plot_de(i), itime, ilev_de(i,:), polyres_de_de)
   ; end do 

   ;  pres=True
   ;  maximize_output(wks2,pres)  ; calls draw and frame for you 


   ntime=dimsizes(intesion)-1
   time=ispan(1,456,1)
 ;****************************************************************************
;随时间变化的量
    diro=diri
    filo1="era_intension.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "定义的各个指数随时间变化的值time=457(1979-201612)"    
    fAtt@creation_date = "2017.8.23"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time"/)
    dimsize=(/ntime/)
    dimunlim=(/False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(time),(/"time"/))
    filevardef(fout, "msf_500", typeof(intesion),(/"time"/))
    filevardef(fout, "STRF", typeof(STRF),(/"time"/))
    filevardef(fout, "slp", typeof(slp),(/"time"/))
    filevardef(fout, "ud", typeof(intension_ud),(/"time"/))
    filevardef(fout, "u10", typeof(u10_mean),(/"time"/))


    filevarattdef(fout, "msf_500", intesion)
    filevarattdef(fout, "STRF", STRF)
    filevarattdef(fout, "slp", slp)
    filevarattdef(fout, "ud", intension_ud)
    filevarattdef(fout, "u10", u10_mean)
    
 
    setfileoption(fout,"DefineMode",False)

 
    fout->time=(/time/)
    fout->msf_500=(/intesion(0:455)/)
    fout->STRF=(/STRF(0:455)/)
    fout->slp=(/slp(0:455)/)
    fout->ud=(/intension_ud(0:455)/)
    fout->u10=(/u10_mean(0:455)/)

end 


