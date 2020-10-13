
;era数据计算的各强度指数和enino3.4的指数的相关性，这里相关性要用enino3.4指数（ano）与各指数的ano做相关。
;计算年循环时候注意时间维的选取，12的整数倍。
;
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
 return intension_an
 end

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


begin
  diri="/mnt/f/share/era/"
  diri1="/mnt/f/share/ersst/"
  fils0=systemfunc("ls "+diri1+"nino3_4.ascii")
  
  nino3_4_r=asciiread(fils0, 457, "float")
  nino3_4=nino3_4_r(0:455)
  delete(nino3_4_r)
  nino3_4!0="time"
  nino3_4&time=ispan(1,456,1)
  printVarSummary(nino3_4)
  siglvl=0.05  ;检验值 
  
  
;-----------------------------------------;
;1.计算定义500hpa以上的msf,经度取140到260强度和位置的相关性 
;-------------------------------------------
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 lon=f1->lon(:)
 msf_time=f1->msf_time
 msf_time_new=where(msf_time.gt.0, msf_time, msf_time@_FillValue)
 copy_VarMeta(msf_time, msf_time_new)
 intension_msf=dim_sum_n_Wrap(msf_time_new(ind(msf_time_new&level.eq.500),ind(lon.eq.140):ind(lon.eq.260),:), 0)  
  intension_msf_ano=clm_ano(intension_msf(0:455))
  r1=escorc_n(intension_msf_ano, nino3_4, 0, 0)
  prob=rtest(r1,dimsizes(intension_msf_ano),0)
 
  yes_no = where(prob.lt.siglvl, True, False)
   print("500hpa以上msf_r1=" + r1 + "yes_no=" +yes_no)
delete(lon)
;*************************************************** 超前滞后相关
 ; mxlag    = 120
 ; x_Lead_y = esccr(intesion,nino3_4,mxlag)
 ; y_Lead_x = esccr(nino3_4,intesion,mxlag)    ; switch the order of the series
 
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
 tind2=ind(date.eq.201612)
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
 slp&time=ispan(1,456,1)
 slp_ano=clm_ano(slp)
 r2=escorc_n(slp_ano, nino3_4, 0, 0)
 prob=rtest(r2,dimsizes(intension_msf_ano),0)
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
 tind2=ind(date.eq.201612)
 
 u10=f3->u10(tind1:tind2,latind1:latind2,lonind1:lonind2)
 u10_mean=dim_avg_n_Wrap(u10, (/1,2/))
 u10_mean_ano=clm_ano(u10_mean)
 r3=escorc_n(u10_mean_ano, nino3_4, 0, 0)

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
 intension_ud_ano =clm_ano(intension_ud(0:455))
 r4=escorc_n(intension_ud_ano, nino3_4, 0, 0)
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
 r5=escorc_n(STRF_ano, nino3_4, 0, 0)
 prob=rtest(r5,dimsizes(STRF_ano),0)
 yes_no = where(prob.lt.siglvl, True, False)
 print("STRF r2=" + r5 + "yes_no=" +yes_no)

;******************************************************************************************
;zero_line location
;**********************************************************************************************
 fils6=systemfunc("ls "+diri+"msf_era_level_location1979-2017.nc")
 f6=addfile(fils6, "r")
 location_xi=f6->locationxi_a(0:455)
 location_xi_ano=clm_ano(location_xi)

 r6=escorc_n(location_xi_ano, nino3_4, 0, 0)
 prob=rtest(r6,dimsizes(location_xi_ano),0)
 yes_no = where(prob.lt.siglvl, True, False)
 print("zero_line r6=" + r6 + "yes_no=" +yes_no)

;******************************************************************************************
;gravity location
;**********************************************************************************************
 fils7=systemfunc("ls "+diri+"msf_era_gravity_center_location1979-2017.nc")
 f7=addfile(fils7, "r")
 x=f7->x(0:455)
 x_ano=clm_ano(x)
 r7=escorc_n(x_ano, nino3_4, 0, 0)
 prob=rtest(r7,dimsizes(x_ano),0)
 yes_no = where(prob.lt.siglvl, True, False)
 print("gravity r7=" + r7 + "yes_no=" +yes_no)
;*************************************************************
 ;超前滞后相关
 ;*******************************************************8
 ; mxlag    = 120
 ; x_Lead_y = esccr(nino3_4,STRF,mxlag)
 ; y_Lead_x = esccr(STRF,nino3_4,mxlag)    ; switch the order of the series
 
 ; ccr = new ( 2*mxlag+1, double)    
 ; ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
 ; ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
 ; time=ispan(-120,120,1)     
 ; wks=gsn_open_wks("png", ""+diri+"位置和强度相关性")
 ; res=True
 ; plot=gsn_csm_xy(wks,time,ccr,res)

;*************************************
 ;各个指数的相关性  intension_msf,slp,u10_mean,intension_ud,STRF
 ;**************************************************************
   ; result=new((/10,2/),double)
   ; result(0,:)=xiangguan(intension_msf,slp)
   ; result(1,:)=xiangguan(intension_msf,u10_mean)
   ; result(2,:)=xiangguan(intension_msf,intension_ud)
   ; result(3,:)=xiangguan(intension_msf,STRF)
   ; result(4,:)=xiangguan(slp,u10_mean)
   ; result(5,:)=xiangguan(slp,intension_ud)
   ; result(6,:)=xiangguan(slp,STRF)
   ; result(7,:)=xiangguan(u10_mean,intension_ud)
   ; result(8,:)=xiangguan(u10_mean,STRF)
   ; result(9,:)=xiangguan(intension_ud,STRF)
   
   ; do i=0,9
   ;  print("result"+i+"=" + result(i,0) +  "yes_no=" +result(i,1))

   ; end do

   ;******************************************************************************
  r_h=(/r2,-r3,r5,r1,r4/)
   r=sprintf("%5.3f",r_h)
   nino3_4_std=runave(dim_standardize_n(nino3_4, 1, 0),13,2)
   slp_ano_std=runave(dim_standardize_n(slp_ano, 1, 0),13,2)
   u10_ano_std=-1*runave(dim_standardize_n(u10_mean_ano, 1, 0),13,2)
   STRF_ano_std=runave(dim_standardize_n(STRF_ano, 1, 0),13,2)
   intension_msf_ano_std=runave(dim_standardize_n(intension_msf_ano, 1, 0),13,2)
   ud_ano_std=runave(dim_standardize_n(intension_ud_ano, 1, 0),13,2)
   location_xi_ano_std=runave(dim_standardize_n(location_xi_ano, 1, 0),13,2)
   x_ano_std=runave(dim_standardize_n(x_ano, 1, 0),13,2)

 
  wks=gsn_open_wks("png", ""+diri+"ersst_nino3_4&指数")
  wks1=gsn_open_wks("png", ""+diri+"ersst_nino3_4&位置指数")

  polyres                  = True     ;;;polyline的属性
  polyres@gsLineThicknessF = 3.0
  polyres@gsLineDashPattern=True


  res = True
  res@gsnFrame = False
  res@gsnDraw = False 

  res@gsnLeftString="(a)dSLP"
  res@gsnRightString="r="+r(0)+""
  res@gsnLeftStringFontHeightF=0.01
  res@gsnRightStringFontHeightF=0.01


  res@gsnAboveYRefLineColor = "gray25"              ; above ref line fill red
  res@gsnBelowYRefLineColor = "gray75"             ; below ref line fill blue
  res@gsnYRefLine           =  0.0             ; create a reference line
  res@tmXBMode="Explicit"
  res@tmXBValues=ispan(1,456,60)
  res@tmXBLabels=ispan(1979,2016,5)
  res@tmXBMinorOn=True
  ; res@tmXBMinorValues=ispan(1,38,1)

  res@tmXBLabelFontAspectF=1.0
  res@tmXBLabelFontHeightF=0.008
  res@tmXBLabelFont=30     ;黑体
  res@tmXBLabelFontThicknessF=2.0
  res@tmXBMajorLengthF=0.005        ;調整刻度的长度
  res@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  ; res@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

  res@tmXBMinorLengthF=0.003
  res@tmXBMinorOutwardLengthF=0.003
 
  res@tmYLMajorLengthF=0.005        ;調整刻度的长度
  res@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  res@tmYLLabelFont=30
  res@tmYLLabelFontHeightF=0.008
  res@tmYLMinorOn=False
  res@tmXTOn=False
  res@tmYROn=False

  res@vpHeightF=0.13
  res@vpWidthF=0.7
  res@vpXF=0.15
  res@vpYF=0.95


plot=new(5,graphic)
lnid=new(5,graphic)
   plot(0)=gsn_csm_xy(wks,nino3_4&time,nino3_4_std,res)
   lnid(0)=gsn_add_polyline(wks,plot(0),nino3_4&time,(/slp_ano_std/),polyres) ; add polyline
    res@gsnLeftString="(b)-U10"
    res@gsnRightString="r="+r(1)+""
    
      res@vpYF=0.77
      
   plot(1)=gsn_csm_xy(wks,nino3_4&time,nino3_4_std,res)
   lnid(1)=gsn_add_polyline(wks,plot(1),nino3_4&time,(/u10_ano_std/),polyres) ; add polyline
    res@gsnLeftString="(c)STRF"
     res@gsnRightString="r="+r(2)+""
 
      res@vpYF=0.59
 
   plot(2)=gsn_csm_xy(wks,nino3_4&time,nino3_4_std,res)
   lnid(2)=gsn_add_polyline(wks,plot(2),nino3_4&time,(/STRF_ano_std/),polyres) ; add polyline

    res@gsnLeftString="(d)msf_500"
     res@gsnRightString="r="+r(3)+""
   
      res@vpYF=0.41
   plot(3)=gsn_csm_xy(wks,nino3_4&time,nino3_4_std,res)
   lnid(3)=gsn_add_polyline(wks,plot(3),nino3_4&time,(/intension_msf_ano_std/),polyres) ; add polyline
   res@gsnLeftString="(e)ud"
    res@gsnRightString="r="+r(4)+""
   
    res@vpYF=0.23

   plot(4)=gsn_csm_xy(wks,nino3_4&time,nino3_4_std,res)
   lnid(4)=gsn_add_polyline(wks,plot(4),nino3_4&time,(/ud_ano_std/),polyres) ; add polyline
   pres=True
    ;  pres@gsnMaximize=True
     
    ; gsn_panel(wks, plot(0:2), (/3,1/), pres)
    maximize_output(wks,pres)  ; calls draw and frame for you 
 
    res@gsnLeftString="(a)msf_zero_line"
    res@gsnRightString="r="+r6+""
    res@vpYF=0.95
    plot_location=new(2,graphic)
    lind_location=new(2,graphic)
    plot_location(0)=gsn_csm_xy(wks1,nino3_4&time,nino3_4_std,res)
    lind_location(0)=gsn_add_polyline(wks1,plot_location(0),nino3_4&time,(/location_xi_ano_std/),polyres) ; add polyline
    res@gsnLeftString="(b)gravity_center"
    res@gsnRightString="r="+r7+""
    res@vpYF=0.77
      
    plot_location(1)=gsn_csm_xy(wks1,nino3_4&time,nino3_4_std,res)
    lind_location(1)=gsn_add_polyline(wks1,plot_location(1),nino3_4&time,(/x_ano_std/),polyres)

    
    pres=True
    pres@gsnMaximize=True
    gsn_panel(wks1, plot_location, (/2,1/), pres) ; calls draw and frame for you 





 end


