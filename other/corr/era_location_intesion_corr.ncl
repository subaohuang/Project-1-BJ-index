;era数据计算的各强度指数和西边界的位置指数的相关性以及各个强度之间的相关性
function   xiangguan(intension,location)
local r1,prob,yes_no,siglvl
begin 
  r1=escorc_n(intension, location, 0, 0)
  prob=rtest(r1,dimsizes(intension),0)
  siglvl=0.05 
  yes_no = where(prob.lt.siglvl, 1.0, 0.0)
   ; print("500hpa以上msf_r1=" + r1 + "yes_no=" +yes_no)
   result=new(2,double)
   result(0)=r1
   result(1)=yes_no
return result
end 

function  an(data1)
begin
 intesion_an=month_to_annual(data1(0:455), 1)
 intesion_an!0="time"
 intesion_an&time=ispan(1,38,1)
 return intesion_an
 end

 function  fitline(data1,data2)
 begin

    rc1=regCoef(data1,data2)
    line = rc1*data1+rc1@yintercept
  return line
  end
;************************************************************
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
;*********************************************************************
begin
 diri="/mnt/f/share/era/"
 fils0=systemfunc("ls "+diri+"msf_era_level_location1979-2017.nc")
 f0=addfile(fils0, "r")
 locationxi_a=f0->locationxi_a
 locationxi_a_ano=clm_ano(locationxi_a(0:455))
 locationdong_a=f0->locationdong_a
 siglvl=0.05  ;检验值
;-------------------------------------------------------;
;1.计算定义500hpa以上的msf,经度取140到260强度和位置的相关性 
;--------------------------------------------------------
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 lon=f1->lon(:)
 msf_time=f1->msf_time
 msf_time_new=where(msf_time.gt.0, msf_time, msf_time@_FillValue)
 copy_VarMeta(msf_time, msf_time_new)
 intension_msf=dim_sum_n_Wrap(msf_time_new(ind(msf_time_new&level.eq.500),ind(lon.eq.140):ind(lon.eq.260),:), 0)  
 intension_msf_ano=clm_ano(intension_msf(0:455))
 r1=escorc_n(intension_msf_ano, locationxi_a_ano, 0, 0)

  prob=rtest(r1,dimsizes(intension_msf_ano),0)
 
  yes_no = where(prob.lt.siglvl, True, False)
   print("500hpa以上msf_r1=" + r1 + "yes_no=" +yes_no)
delete(lon)
;*************************************************** 超前滞后相关
 ; mxlag    = 120
 ; x_Lead_y = esccr(intesion,locationxi_a,mxlag)
 ; y_Lead_x = esccr(locationxi_a,intesion,mxlag)    ; switch the order of the series
 
 ; ccr = new ( 2*mxlag+1, double)    
 ; ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
 ; ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
 ; x=ispan(-120,120,1)     
 ; wks=gsn_open_wks("png", ""+diri+"位置和强度相关性")
 ; res=True
 ; plot=gsn_csm_xy(wks,x,ccr,res)
 ; delete(lon)
 ; delete(level)
 ; delete(f1)
   
;*******************************************************************
;-----------------------------------------------------------------
;2.slp 与位置相关
;-----------------------------------------------------------------

 fils2=systemfunc("ls "+diri+"era_slp_80E-70w.nc")
 f2=addfile(fils2, "r") 

 time=f2->time(:)
 lon=f2->longitude(:)
 lat=f2->latitude(:)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 date=cd_calendar(time, 1)
 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201701)
 msl=f2->msl(tind1:tind2,latind1:latind2,:)
 m=getvaratts(msl)

  if (any(getvaratts(msl).eq."scale_factor")) then 
      msl_new=dble2flt(msl*msl@scale_factor+msl@add_offset)
  end if
  copy_VarCoords(msl, msl_new)
  delete(msl)
 msl_tahiti=dim_avg_n_Wrap(msl_new(latitude|:,longitude|ind(lon.eq.200):ind(lon.eq.280),time|:),(/0,1/))
 msl_daerwen=dim_avg_n_Wrap(msl_new(latitude|:,longitude|ind(lon.eq.80):ind(lon.eq.160),time|:),(/0,1/))
 slp=msl_tahiti - msl_daerwen

 slp!0="time"
 slp&time=ispan(1,457,1)
 slp_ano=clm_ano(slp(0:455))
 r2=escorc_n(slp_ano, locationxi_a_ano, 0, 0)
 prob=rtest(r2,dimsizes(intension_msf),0)
  yes_no = where(prob.lt.siglvl, True, False)
 print("deta slp r2=" + r2 + "yes_no=" +yes_no)
  
delete(lon)
delete(lat)
delete(time)
delete(date)
;******************************************************************************************
;us-1和位置的关系 5度到-5度，150E到150W的结果。
;**********************************************************************************************
 fils3=systemfunc("ls "+diri+"era_u_sst_slp.nc")
 f3=addfile(fils3, "r")

 time=f3->time(:)
 lon=f3->longitude(:)
 lat=f3->latitude(:)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 lonind1=ind(lon.eq.150)
 lonind2=ind(lon.eq.210)
 
 date=cd_calendar(time, 1)
 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201701)
 
 u10=f3->u10(tind1:tind2,latind1:latind2,lonind1:lonind2)
  if (any(getvaratts(u10).eq."scale_factor")) then 
      u10_new=dble2flt(u10*u10@scale_factor+u10@add_offset)
  end if
  copy_VarCoords(u10, u10_new)
  delete(u10)
 u10_mean=dim_avg_n_Wrap(u10_new, (/1,2/))
 u10_mean_ano=clm_ano(u10_mean(0:455))
 r3=escorc_n(u10_mean_ano, locationxi_a_ano, 0, 0)

 prob=rtest(r3,dimsizes(intension_msf_ano),0)
  yes_no = where(prob.lt.siglvl, True, False)

  print("u10 r3=" + r3 +"yes_no=" +yes_no)
delete(time)
delete(lon)
delete(lat)
delete(date)
;******************************************************************************************
;ud 和位置的关系,ud定义的强度，基于平均的结果，150E到210度E，150到300和1000到750的平均结果
;**********************************************************************************************
 fils4=systemfunc("ls "+diri+"msf_time_era.nc")
 f4=addfile(fils4, "r")
 level=f4->level
 lon=f4->lon
 lonind1=ind(lon.eq.150)
 lonind2=ind(lon.eq.210)
 ud_time=f4->ud_time 
 intension_ud=dim_avg_n_Wrap(ud_time(ind(level.eq.300):ind(level.eq.150),lonind1:lonind2,:),(/0,1/))-dim_avg_n_Wrap(ud_time(ind(level.eq.1000):ind(level.eq.750),lonind1:lonind2,:),(/0,1/))

 intension_ud!0="time"
 intension_ud&time=ispan(1,457,1)
 intension_ud_ano=clm_ano(intension_ud(0:455))
 r4=escorc_n(intension_ud_ano, locationxi_a_ano, 0, 0)
 prob=rtest(r4,dimsizes(intension_ud_ano),0)
  yes_no = where(prob.lt.siglvl, True, False)
 print("ud r2=" + r4 + "yes_no=" +yes_no)
 delete(level)
delete(lon)
;******************************************************************************************
;STRF
;**********************************************************************************************
 fils5=systemfunc("ls "+diri+"msf_time_era.nc")
 f5=addfile(fils5, "r")
 level=f5->level(:)
 lon=f5->lon(:)
 lonind1=ind(lon.eq.150)
 lonind2=ind(lon.eq.210)
 msf_time=f5->msf_time
 STRF=dim_avg_n_Wrap(msf_time(:,lonind1:lonind2,:),(/0,1/))
 STRF_ano=clm_ano(STRF(0:455))
 r5=escorc_n(STRF_ano, locationxi_a_ano, 0, 0)
 prob=rtest(r5,dimsizes(STRF_ano),0)
 yes_no = where(prob.lt.siglvl, True, False)
 print("STRF r2=" + r5 + "yes_no=" +yes_no)
 ;*****************************************************************
 ;各个指数的相关性  intension_msf,slp,u10_mean,intension_ud,STRF
 ;**************************************************************
   result=new((/10,2/),double)
   result(0,:)=xiangguan(intension_msf_ano,slp_ano)
   result(1,:)=xiangguan(intension_msf_ano,u10_mean_ano)
   result(2,:)=xiangguan(intension_msf_ano,intension_ud_ano)
   result(3,:)=xiangguan(intension_msf_ano,STRF_ano)
   result(4,:)=xiangguan(slp_ano,u10_mean_ano)
   result(5,:)=xiangguan(slp_ano,intension_ud_ano)
   result(6,:)=xiangguan(slp_ano,STRF_ano)
   result(7,:)=xiangguan(u10_mean_ano,intension_ud_ano)
   result(8,:)=xiangguan(u10_mean_ano,STRF_ano)
   result(9,:)=xiangguan(intension_ud_ano,STRF_ano)
   
   do i=0,9
    print("result"+i+"=" + result(i,0) +  "yes_no=" +result(i,1))

   end do
  ;*********************************************************************
  ;画强度与位置的散点图，表明情况
  ;************************************************************************
    locationxi_a_an=an(locationxi_a)
    intension_slp_an=an(slp)

    intension_u10_an=-1*an(u10_mean)
    copy_VarCoords(intension_slp_an,intension_u10_an)

    intension_STRF_an=an(STRF)

    intension_msf_an=an(intension_msf)
    
    intension_ud_an=an(intension_ud)

    line1=fitline(locationxi_a_an,intension_slp_an)
    line2=fitline(locationxi_a_an,intension_u10_an)
    line3=fitline(locationxi_a_an,intension_STRF_an)
    line4=fitline(locationxi_a_an,intension_msf_an)
    line5=fitline(locationxi_a_an,intension_ud_an)
    
    

   wks=gsn_open_wks("png", ""+diri+"scatter_intension&location")

   resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@xyLineThicknesses = (/ 1.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"Blue"  ,"black"/)          ; change line color
    resxy@xyMarkerColor = "Blue"
    resxy@xyMarkLineModes   = (/"Markers","lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    resxy@xyMarkers         =  16                      ; choose type of marker  
    resxy@xyMarkerSizeF     = 0.005        ; Marker size (default 0.01)
    
    
    resxy@gsnLeftString="(a)dSLP"
    resxy@gsnRightString="msf_zeroline"
    ; resxy@gsnRightString=""+slp_data(3,1)+""
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01
    ; resxy@trXMaxF=dimsizes(intesion&time)
    ; resxy@tmXBValues=ispan(13,dimsizes(intesion&time),60)
    ; resxy@tmXBLabels=ispan(1980,2015,5)

    resxy@vpHeightF=0.13
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.95
    
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
    resxy@tmYLMinorOn=False
    resxy@tmXTOn=False
    resxy@tmYROn=False
    

   plot=new(5,graphic)

    
   plot(0)=gsn_csm_xy(wks,locationxi_a_an,(/intension_slp_an,line1/),resxy)
    resxy@gsnLeftString="(b)-U10"
    ; resxy@gsnRightString=""+u10_data(3,1)+""
      resxy@vpYF=0.77
     
   plot(1)=gsn_csm_xy(wks,locationxi_a_an,(/intension_u10_an,line2/),resxy)
    resxy@gsnLeftString="(c)STRF"
    ; resxy@gsnRightString=""+STRF_data(3,1)+""
      resxy@vpYF=0.59
      resxy@tmYLMode="Explicit"    
       resxy@tmYLValues=fspan(0.4*10^9, 3.2*10^9, 8)
      resxy@tmYLLabels=fspan(0.4,3.2,8)

  plot(2)=gsn_csm_xy(wks,locationxi_a_an,(/intension_STRF_an,line3/),resxy)
      delete(resxy@tmYLLabels)
      delete(resxy@tmYLValues)
      resxy@gsnLeftString="(d)msf_500"
    ; resxy@gsnRightString=""+data_an(3,1)+""
      resxy@vpYF=0.41
       resxy@tmYLValues=fspan(0.8*10^11, 2.0*10^11, 7)
       resxy@tmYLLabels=fspan(0.8,2.0,7)
  plot(3)=gsn_csm_xy(wks,locationxi_a_an,(/intension_msf_an,line4/),resxy)
     resxy@gsnLeftString="(e)ud"
   ; resxy@gsnRightString=""+ud_data(3,1)+""
     resxy@vpYF=0.23
     delete(resxy@tmYLMode)
     delete(resxy@tmYLLabels)
     delete(resxy@tmYLValues)
   plot(4)=gsn_csm_xy(wks,locationxi_a_an,(/intension_ud_an,line5/),resxy)
      
   

txid_tr = new(5, graphic)
amid_tr = new(5, graphic)
txres                       = True
txres@txPerimOn             = True
txres@txFontHeightF         = 0.01

amres_tr                  = True
amres_tr@amParallelPosF   =  0.45   ; This is the right edge of the plot.
amres_tr@amOrthogonalPosF =  -0.35    ; This is the top edge of the plot.

r=(/r2,r3,r5,r1,r4/)
letter = sprintf("%5.3f",r)
do i = 0, 4
;---Create text strings
    tr_label = "("  +"r="+ letter(i) + ")" 

    txres@txBackgroundFillColor = "transparent"
    txres@txPerimOn = False
    txid_tr(i) = gsn_create_text(wks, tr_label, txres)

;---Attach text strings to plot
    amid_tr(i) = gsn_add_annotation(plot(i), txid_tr(i), amres_tr)

end do

  pres=True
    ;  pres@gsnMaximize=True
    ; gsn_panel(wks, plot(0:2), (/3,1/), pres)
    maximize_output(wks,pres)  ; calls draw and frame for you 
    
 end
