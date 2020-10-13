
 ;******************************************************************
 ;enino modoki msf
;***************************************************************
function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin
    ; data_new=doubletofloat(dim_standardize_n_Wrap(data(::-1,:,:),1,(/0,1/)))       ;标准化
    ; print(data_new)
    ; copy_VarMeta(data(::-1,:,:), data_new)
    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,2)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
     ; print(rc)
return(rc)
end
;--------------------------------------------------------------------------
function modoki_season_qcx(modoki_yr,season)  ;chun :2-4,取2，夏：5 秋:8 冬：10,modoki_yr指的是秋季发生的年份的下一年 
                                                         ;
begin
month=ispan(1,12,1)
 if season.ne.11 then
 modoki_month_sp=new((/dimsizes(modoki_yr),3/),string)
  do i=0,dimsizes(modoki_yr)-1
    do j=season,season+2
  modoki_month_sp(i,j-season)=""+modoki_yr(i)+""+""+month(j)+""
  end do
  end do
 end if
modoki_month_sp_new=ndtooned ( stringtointeger(modoki_month_sp))
return modoki_month_sp_new
end 
;------------------------------------------------------
function modoki_season_dong(modoki_yr1,modoki_yr2)   ;chun :2-4,取2，夏：5 秋:8 冬：10,modoki_yr指的是秋季发生的年份的下一年 
                                                         ;
begin
month=ispan(1,12,1)
 
 modoki_month_sp=new((/dimsizes(modoki_yr1),3/),string)
  do i=0,dimsizes(modoki_yr1)-1
    do j=1,2
  modoki_month_sp(i,j)=""+modoki_yr2(i)+""+""+month(j-1)+""
 end do
 modoki_month_sp(i,0)=""+modoki_yr1(i)+""+""+month(11)+""
end do
modoki_month_sp_new=ndtooned ( stringtointeger(modoki_month_sp))
return modoki_month_sp_new
end 
;----------------------------------------------------------
function modoki_anu(modoki_yr1,modoki_yr2)
begin
 month=ispan(1,12,1)
 modoki_month_sp=new((/dimsizes(modoki_yr1),12/),string)
  do i=0,dimsizes(modoki_yr1)-1
    do j=0,9
  modoki_month_sp(i,j)=""+modoki_yr1(i)+""+""+month(j+2)+""
  end do
    do j=10,11
  modoki_month_sp(i,j)=""+modoki_yr2(i)+""+""+month(j-10)+""
  end do
end do
modoki_month_sp_new=ndtooned ( stringtointeger(modoki_month_sp))
return modoki_month_sp_new
end 

;-------------------------------------------------------------
function tiqu_modoki(msf_scale, ud_time, w_time,  yr_month , modoki_yrI , modoki_yrI1)
begin

modokiI_month_sp=modoki_season_qcx(modoki_yrI,2)
modokiI_month_sum=modoki_season_qcx(modoki_yrI,5)
modokiI_month_atu=modoki_season_qcx(modoki_yrI,8)
modokiI_month_win=modoki_season_dong(modoki_yrI,modoki_yrI1)

ind_modokiI_month_sp=get1Dindex(yr_month,modokiI_month_sp)
ind_modokiI_month_sum=get1Dindex(yr_month,modokiI_month_sum)
ind_modokiI_month_atu=get1Dindex(yr_month,modokiI_month_atu)
ind_modokiI_month_win=get1Dindex(yr_month,modokiI_month_win)

msf_modokiI_sp=msf_scale(:,:,ind_modokiI_month_sp)
msf_modokiI_sum=msf_scale(:,:,ind_modokiI_month_sum)
msf_modokiI_atu=msf_scale(:,:,ind_modokiI_month_atu)
msf_modokiI_win=msf_scale(:,:,ind_modokiI_month_win)

msf_modokiI_sp_avg=dim_avg_n_Wrap( msf_modokiI_sp,2)
msf_modokiI_sum_avg=dim_avg_n_Wrap( msf_modokiI_sum,2)
msf_modokiI_atu_avg=dim_avg_n_Wrap( msf_modokiI_atu,2)
msf_modokiI_win_avg=dim_avg_n_Wrap( msf_modokiI_win,2)


ud_modokiI_sp_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modokiI_month_sp),2)
ud_modokiI_sum_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modokiI_month_sum),2)
ud_modokiI_atu_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modokiI_month_atu),2)
ud_modokiI_win_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modokiI_month_win),2)

w_modokiI_sp_avg=dim_avg_n_Wrap( w_time(:,:,ind_modokiI_month_sp),2)
w_modokiI_sum_avg=dim_avg_n_Wrap( w_time(:,:,ind_modokiI_month_sum),2)
w_modokiI_atu_avg=dim_avg_n_Wrap( w_time(:,:,ind_modokiI_month_atu),2)
w_modokiI_win_avg=dim_avg_n_Wrap( w_time(:,:,ind_modokiI_month_win),2)
   
msf_modokiI=new((/4,dimsizes(msf_scale&level),dimsizes(msf_scale&lon)/),typeof(msf_scale))
ud_modokiI=new((/4,dimsizes(msf_scale&level),dimsizes(msf_scale&lon)/),typeof(msf_scale))
w_modokiI=new((/4,dimsizes(msf_scale&level),dimsizes(msf_scale&lon)/),typeof(msf_scale))
modokiI=new(4,float)

msf_modokiI(0,:,:)=msf_modokiI_sp_avg
msf_modokiI(1,:,:)=msf_modokiI_sum_avg
msf_modokiI(2,:,:)=msf_modokiI_atu_avg
msf_modokiI(3,:,:)=msf_modokiI_win_avg

ud_modokiI(0,:,:)=ud_modokiI_sp_avg
ud_modokiI(1,:,:)=ud_modokiI_sum_avg
ud_modokiI(2,:,:)=ud_modokiI_atu_avg
ud_modokiI(3,:,:)=ud_modokiI_win_avg


w_modokiI(0,:,:)=w_modokiI_sp_avg
w_modokiI(1,:,:)=w_modokiI_sum_avg
w_modokiI(2,:,:)=w_modokiI_atu_avg
w_modokiI(3,:,:)=w_modokiI_win_avg


modokiI@msf_modokiI=msf_modokiI
modokiI@ud_modokiI=ud_modokiI
modokiI@w_modokiI=w_modokiI

return modokiI
end 

;****************************************************************
begin
 diri="/mnt/f/share/era/data/"
 diri_pi="/mnt/f/share/era/"
 fils1=systemfunc("ls "+diri+"msf_time_era_global_10--10.nc")
 f1=addfile(fils1, "r")
 msf_time=f1->msf_time(:,:,0:455)
 ud_time_r=f1->ud_time(:,:,0:455)
 w_time_r=f1->w_time(:,:,0:455)
 level=f1->level(:)
 lon=f1->lon(:)
 
 msf_time&level@units="hpa"
 msf_scale=dtrend_n(rmMonAnnCycLLT(msf_time),True,2)/10^9
 copy_VarMeta(msf_time, msf_scale)
 ud_time=dtrend_n(rmMonAnnCycLLT(ud_time_r),True,2)
 w_time=dtrend_n(rmMonAnnCycLLT(w_time_r),True,2)*5
  copy_VarMeta(ud_time_r, ud_time)
  copy_VarMeta(w_time_r, w_time)
 ;---------------------------------------------------------------
 month=ispan(1,12,1)
 yr=ispan(1979,2016,1) 
 yr_month=new(dimsizes(msf_time&time),string)
do t=0,455
if t.le.11 then 
h=mod(t,12)
    yr_month(t)=""+yr(0)+""+""+month(t)+""
end if   

if t.gt.11 then 
h=mod(t,12)
m=(t-h)/12
    yr_month(t)=""+yr(m)+""+""+month(h)+""
end if
end do 
;----------------------------------------------

;------------------------------------------------------------------------
;挑选enino modoki的月份，3,4,5春季，6,7,8为夏季，9,10,11为秋季，12,1,2为冬季.
;------------------------------------------------------------------------
   
 modoki_yrI=(/1987,1990,2002/)
 modoki_yrII=(/1979,1991,1992,2004,2009/)
 canonical_yr=(/1982,1997/)

 modoki_yrI1=modoki_yrI+1
 modoki_yrII1=modoki_yrII+1
 canonical_yr1=canonical_yr+1
 

;---------------------------------------------------------------------------

   
 modokiI=tiqu_modoki(msf_scale, ud_time, w_time,  yr_month , modoki_yrI , modoki_yrI1)
 msf_modokiI=modokiI@msf_modokiI
 ud_modokiI=modokiI@ud_modokiI
 w_modokiI=modokiI@w_modokiI

 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), msf_modokiI)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), ud_modokiI)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), w_modokiI)

 msf_modokiI_sp=msf_modokiI(0,:,:)
 msf_modokiI_sum=msf_modokiI(1,:,:)
 msf_modokiI_atu=msf_modokiI(2,:,:)
 msf_modokiI_win=msf_modokiI(3,:,:)


 ud_modokiI_sp=ud_modokiI(0,:,:)
 ud_modokiI_sum=ud_modokiI(1,:,:)
 ud_modokiI_atu=ud_modokiI(2,:,:)
 ud_modokiI_win=ud_modokiI(3,:,:)

 w_modokiI_sp=w_modokiI(0,:,:)
 w_modokiI_sum=w_modokiI(1,:,:)
 w_modokiI_atu=w_modokiI(2,:,:)
 w_modokiI_win=w_modokiI(3,:,:)


   
 modokiII=tiqu_modoki(msf_scale, ud_time, w_time,  yr_month , modoki_yrII , modoki_yrII1)
 msf_modokiII=modokiII@msf_modokiI
 ud_modokiII=modokiII@ud_modokiI
 w_modokiII=modokiII@w_modokiI

 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), msf_modokiII)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), ud_modokiII)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), w_modokiII)

 msf_modokiII_sp=msf_modokiII(0,:,:)
 msf_modokiII_sum=msf_modokiII(1,:,:)
 msf_modokiII_atu=msf_modokiII(2,:,:)
 msf_modokiII_win=msf_modokiII(3,:,:)


 ud_modokiII_sp=ud_modokiII(0,:,:)
 ud_modokiII_sum=ud_modokiII(1,:,:)
 ud_modokiII_atu=ud_modokiII(2,:,:)
 ud_modokiII_win=ud_modokiII(3,:,:)

 w_modokiII_sp=w_modokiII(0,:,:)
 w_modokiII_sum=w_modokiII(1,:,:)
 w_modokiII_atu=w_modokiII(2,:,:)
 w_modokiII_win=w_modokiII(3,:,:)


 canonical=tiqu_modoki(msf_scale, ud_time, w_time,  yr_month , canonical_yr , canonical_yr1)
 msf_canonical = canonical@msf_modokiI
 ud_canonical  = canonical@ud_modokiI
 w_canonical = canonical@w_modokiI

 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), msf_canonical)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), ud_canonical)
 copy_VarCoords(msf_scale(time|0:3,level|:,lon|:), w_canonical)

 msf_canonical_sp=msf_canonical(0,:,:)
 msf_canonical_sum=msf_canonical(1,:,:)
 msf_canonical_atu=msf_canonical(2,:,:)
 msf_canonical_win=msf_canonical(3,:,:)


 ud_canonical_sp=ud_canonical(0,:,:)
 ud_canonical_sum=ud_canonical(1,:,:)
 ud_canonical_atu=ud_canonical(2,:,:)
 ud_canonical_win=ud_canonical(3,:,:)


 w_canonical_sp=w_canonical(0,:,:)
 w_canonical_sum=w_canonical(1,:,:)
 w_canonical_atu=w_canonical(2,:,:)
 w_canonical_win=w_canonical(3,:,:)
 


 modoki1_anu = modoki_anu(modoki_yrI,modoki_yrI1)
 modoki2_anu = modoki_anu(modoki_yrII,modoki_yrII1)
 canonical_anu=modoki_anu(canonical_yr,canonical_yr1)
 
 ind_modokiI_anu=get1Dindex(yr_month,modoki1_anu)
 ind_modokiII_anu=get1Dindex(yr_month,modoki2_anu)
 ind_canonical_anu=get1Dindex(yr_month,canonical_anu)


 msf_modokiI_anu=dim_avg_n_Wrap(msf_scale(:,:,ind_modokiI_anu),2)
 msf_modokiII_anu=dim_avg_n_Wrap(msf_scale(:,:,ind_modokiII_anu),2)
 msf_canonical_anu=dim_avg_n_Wrap(msf_scale(:,:,ind_canonical_anu),2)

 ud_modokiI_anu=dim_avg_n_Wrap(ud_time(:,:,ind_modokiI_anu),2)
 ud_modokiII_anu=dim_avg_n_Wrap(ud_time(:,:,ind_modokiII_anu),2)
 ud_canonical_anu=dim_avg_n_Wrap(ud_time(:,:,ind_canonical_anu),2)

 w_modokiI_anu=dim_avg_n_Wrap(w_time(:,:,ind_modokiI_anu),2) 
 w_modokiII_anu=dim_avg_n_Wrap(w_time(:,:,ind_modokiII_anu),2) 
 w_canonical_anu=dim_avg_n_Wrap(w_time(:,:,ind_canonical_anu),2) 

printVarSummary(w_modokiI_anu)


;**********************************************************************
 ;画图
;************************************************************************
  wks=gsn_open_wks("pdf", ""+diri_pi+"era_modoki异常值的四季_10到-10平均的结果")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"

    
    ; res@cnFillPalette  = "BlueDarkRed18"
    ; gsn_define_colormap(wks,"hotcolr_19lev")
    gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -2
    res@cnMaxLevelValF = 2
    res@cnLevelSpacingF = 0.5
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    ; res@cnLinePalette="BlueDar0kRed18"
    ; res@cnLinePalette="hotcolr_19lev"
    res@cnLineLabelsOn=True
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False  
    res@tiMainString="canonical El Nino"
    ; res@cnFillColors=(/2,4,5,6,10,10,14,16,17,20/)
     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
     ; res@cnLabelScalingMode="ScaleFactor"
     ; res@cnLabelScaleValueF=10^9
      
    res@lbLabelBarOn=False
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
      ; res@tmXBMode="Explicit"
      ; res@tmXBValues=ispan(100,290,40)
      ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
      ; res@tmXBMinorOn=False

    res@tmXBMode="Explicit"
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)



    res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    res@vcRefLengthF    = 0.04        ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.02       ; thin out vectors
    res@vcMinDistanceF  = 0.03       ; thin out vectors
    res@vcMapDirection  = False

    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15
    
    res@vcRefAnnoOrthogonalPosF  = -1.05
    res@vcRefAnnoPerimOn=False
    ; res@vcGlyphStyle="CurlyVector"
    res@vcGlyphStyle="FillArrow"
    res@vcLineArrowColor="black"
    res@vcLineArrowThicknessF=1.5
    res@tmYRMode ="Automatic"

    plot=new(12,graphic)

   plot(0)=gsn_csm_pres_hgt_vector(wks,msf_canonical_sp,ud_canonical_sp,w_canonical_sp,res)
    res@tiMainString="El Nino Modoki I"
     res@vcRefMagnitudeF = 1.5              
   plot(1)=gsn_csm_pres_hgt_vector(wks,msf_modokiI_sp,ud_modokiI_sp,w_modokiI_sp,res)
    res@tiMainString="El Nino Modoki II"
     res@vcRefMagnitudeF = 1.5              
   plot(2)=gsn_csm_pres_hgt_vector(wks,msf_modokiII_sp,ud_modokiII_sp,w_modokiII_sp,res)
    delete(res@tiMainString)
    res@vcRefMagnitudeF = 2.0 
   plot(3)=gsn_csm_pres_hgt_vector(wks,msf_canonical_sum,ud_canonical_sum,w_canonical_sum,res)
    res@vcRefMagnitudeF = 1.5
   plot(4)=gsn_csm_pres_hgt_vector(wks,msf_modokiI_sum,ud_modokiI_sum,w_modokiI_sum,res)
    res@vcRefMagnitudeF = 1.5
   plot(5)=gsn_csm_pres_hgt_vector(wks,msf_modokiII_sum,ud_modokiII_sum,w_modokiII_sum,res)
     res@vcRefMagnitudeF = 2.0 
   plot(6)=gsn_csm_pres_hgt_vector(wks,msf_canonical_atu,ud_canonical_atu,w_canonical_atu,res)
     res@vcRefMagnitudeF = 1.5
   plot(7)=gsn_csm_pres_hgt_vector(wks,msf_modokiI_atu,ud_modokiI_atu,w_modokiI_atu,res)
    res@vcRefMagnitudeF = 1.5
   plot(8)=gsn_csm_pres_hgt_vector(wks,msf_modokiII_atu,ud_modokiII_atu,w_modokiII_atu,res)
    res@vcRefMagnitudeF = 2.0 
   plot(9)=gsn_csm_pres_hgt_vector(wks,msf_canonical_win,ud_canonical_win,w_canonical_win,res)
    res@vcRefMagnitudeF = 1.5
   plot(10)=gsn_csm_pres_hgt_vector(wks,msf_modokiI_win,ud_modokiI_win,w_modokiI_win,res)
     res@vcRefMagnitudeF = 1.5
   plot(11)=gsn_csm_pres_hgt_vector(wks,msf_modokiII_win,ud_modokiII_win,w_modokiII_win,res)
                                                 
;-------------------------------
delete(res@gsnRightString)
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 3.
res@cnFillOn        = False
res@tmYRMode ="Automatic"

plot1=new(12,graphic)
   plot1(0)=gsn_csm_pres_hgt(wks,msf_canonical_sp,res)
    
   plot1(1)=gsn_csm_pres_hgt(wks,msf_modokiI_sp,res)
  
   plot1(2)=gsn_csm_pres_hgt(wks,msf_modokiII_sp,res)
  
   plot1(3)=gsn_csm_pres_hgt(wks,msf_canonical_sum,res)
    
   plot1(4)=gsn_csm_pres_hgt(wks,msf_modokiI_sum,res)
  
   plot1(5)=gsn_csm_pres_hgt(wks,msf_modokiII_sum,res)

   plot1(6)=gsn_csm_pres_hgt(wks,msf_canonical_atu,res)
    
   plot1(7)=gsn_csm_pres_hgt(wks,msf_modokiI_atu,res)
  
   plot1(8)=gsn_csm_pres_hgt(wks,msf_modokiII_atu,res)

   plot1(9)=gsn_csm_pres_hgt(wks,msf_canonical_win,res)
    
   plot1(10)=gsn_csm_pres_hgt(wks,msf_modokiI_win,res)
  
   plot1(11)=gsn_csm_pres_hgt(wks,msf_modokiII_win,res)

;---Resources for the polyline  添加下垫面的陆地
  ployres                  = True          ; polyline mods desired
  ployres@gsLineThicknessF = 15.0            ; line thickness   
  ployres@gsLineColor      = "black"         ; color of lines
  ployres@gsLineDashPattern=0
  ;---Attach the polyline
  ploylevel=1000
  ploylon1=fspan(9, 39, 25)
  ploylon2=fspan(100, 149, 25)
  ploylon3=fspan(280, 315, 25)
  ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
  ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
  ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)

  dum1=new(12,graphic)
  dum2=new(12,graphic) 
  dum3=new(12,graphic) 
  do i=0,11
 dum1(i)= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
  end do
ployres@gsLineDashPattern=15
do i=0,11
 dum2(i)= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel1,  ployres)
  end do
  ployres@gsLineDashPattern=0
do i=0,11
 dum3(i)= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel1,  ployres)
  end do


do i=0,11
overlay(plot(i),plot1(i))
end do 
pres=True
pres@gsnPanelLabelBar=True
pres@gsnMaximize=True
pres@gsnPanelFigureStrings=(/"MMA","MMA","MMA","JJA","JJA","JJA","SON","SON","SON","DJF","DJF","DJF"/) 
pres@gsnPanelFigureStringsJust="TopLeft"
gsn_panel(wks, plot, (/4,3/), pres)
delete(pres)
;---------------------------------------------------------
;画年平均的变化
;--------------------------------------------------------
wks2=gsn_open_wks("png", ""+diri_pi+"era_modoki异常值的年平均_10--10")
delete(res@cnLineThicknessF)
delete(res@vcRefMagnitudeF)
delete(res@cnLevels)

gsn_define_colormap(wks2,"hotcolr_19lev")
delete(res@cnFillColors)
    res@cnFillOn        = True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -3
    res@cnMaxLevelValF = 3
    res@cnLevelSpacingF = 0.5
    res@cnLineLabelsOn=False
   res@cnFillColors=(/2,3,4,5,6,7,10,10,14,15,16,17,20,21/)
   plot_anu=new(3,graphic)
   res@tiMainString="canonical El Nino"
   plot_anu(0)=gsn_csm_pres_hgt_vector(wks2,msf_canonical_anu,ud_canonical_anu,w_canonical_anu,res)
   delete(res@cnFillColors)
    res@cnMinLevelValF = -2
    res@cnMaxLevelValF = 2
    res@cnLevelSpacingF = 0.5
    res@cnFillColors=(/2,4,5,6,10,10,14,16,17,20/)
    res@tiMainString="El Nino Modoki I"
     res@vcRefMagnitudeF = 1.0              
   plot_anu(1)=gsn_csm_pres_hgt_vector(wks2,msf_modokiI_anu,ud_modokiI_anu,w_modokiI_anu,res)
    res@tiMainString="El Nino Modoki II"
     res@vcRefMagnitudeF = 1.0              
   plot_anu(2)=gsn_csm_pres_hgt_vector(wks2,msf_modokiII_anu,ud_modokiII_anu,w_modokiII_anu,res)
    delete(res@tiMainString)
    res@vcRefMagnitudeF = 2.0 
    
  ;-------------------------------
 
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 3.
res@cnFillOn        = False
res@tmYRMode ="Automatic"

plot1_anu=new(3,graphic)
   plot1_anu(0)=gsn_csm_pres_hgt(wks2,msf_canonical_anu,res)
    
   plot1_anu(1)=gsn_csm_pres_hgt(wks2,msf_modokiI_anu,res)
  
   plot1_anu(2)=gsn_csm_pres_hgt(wks2,msf_modokiII_anu,res)

  dum4=new(3,graphic)
  dum5=new(3,graphic) 
  dum6=new(3,graphic)  
  do i=0,2
 dum4(i)= gsn_add_polyline(wks2,plot_anu(i),  ploylon1,ploylevel1,  ployres)
  end do
ployres@gsLineDashPattern=15
do i=0,2
 dum5(i)= gsn_add_polyline(wks2,plot_anu(i),  ploylon2,ploylevel1,  ployres)
  end do
 ployres@gsLineDashPattern=0
do i=0,2
 dum6(i)= gsn_add_polyline(wks2,plot_anu(i),  ploylon3,ploylevel1,  ployres)
  end do
;;;

do i=0,2
overlay(plot_anu(i),plot1_anu(i))
end do 
pres=True
; pres@gsnPanelLabelBar=True
pres@gsnMaximize=True
gsn_panel(wks2, plot_anu, (/1,3/), pres)
; k=plot_anu(0)
; draw(k)
; frame(wks2)
end
