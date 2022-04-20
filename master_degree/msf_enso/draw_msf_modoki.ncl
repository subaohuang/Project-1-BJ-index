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

;****************************************************************
begin
 diri="/mnt/f/share/era/"
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 msf_time=f1->msf_time(:,:,:)
 ud_time=f1->ud_time(:,:,:)
 w_time=f1->w_time(:,:,:)
 level=f1->level(:)
 lon=f1->lon(:)
 msf_time&level@units="hpa"
 msf_scale=msf_time/10^9
 copy_VarMeta(msf_time, msf_scale)

 printVarSummary(ud_time)
;------------------------------------------------------------------------
;挑选enino modoki的月份，3,4,5春季，6,7,8为夏季，9,10,11为秋季，12,1,2为冬季.
;------------------------------------------------------------------------
   
 modoki_yrI=(/1979,1987,1990,1991,1992,1994,2002,2004,2009/)
 modoki_yr2=modoki_yr+1

 month=ispan(1,12,1)
 yr=ispan(1979,2017,1) 
 yr_month=new(dimsizes(msf_time&time),string)


;---------------------------------------------------------------------------
;转化时间
;  转化时间维197901--------------------------------------------------------------------------
 ; do t=0,456
 ;    if t.le.11 then 
 ;    h=mod(t,12)
 ;      if t.le.8 then
 ;        yr_month(t)=""+yr(h)+""+"0"+month(t)+""
 ;      else
 ;       yr_month(t)=""+yr(h)+""+""+month(t)+""
 ;  end if  
 ;  end if   

 ;    if t.gt.11 then 
 ;   h=mod(t,12)
 ;   m=(t-h)/12
 ;   if h.le.8 then 
 ;        yr_month(t)=""+yr(m)+""+"0"+month(h)+""
 ;  else 
 ;        yr_month(t)=""+yr(m)+""+""+month(h)+""
 ;    end if
 ;   end if
 ;  end do 
 ;-----------------------------------------------------------------------------

do t=0,456
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

modoki_month_sp=modoki_season_qcx(modoki_yr,2)
modoki_month_sum=modoki_season_qcx(modoki_yr,5)
modoki_month_atu=modoki_season_qcx(modoki_yr,8)
modoki_month_win=modoki_season_dong(modoki_yr,modoki_yr2)

ind_modoki_month_sp=get1Dindex(yr_month,modoki_month_sp)
ind_modoki_month_sum=get1Dindex(yr_month,modoki_month_sum)
ind_modoki_month_atu=get1Dindex(yr_month,modoki_month_atu)
ind_modoki_month_win=get1Dindex(yr_month,modoki_month_win)

msf_modoki_sp=msf_scale(:,:,ind_modoki_month_sp)
msf_modoki_sum=msf_scale(:,:,ind_modoki_month_sum)
msf_modoki_atu=msf_scale(:,:,ind_modoki_month_atu)
msf_modoki_win=msf_scale(:,:,ind_modoki_month_win)

msf_modoki_sp_avg=dim_avg_n_Wrap( msf_modoki_sp,2)
msf_modoki_sum_avg=dim_avg_n_Wrap( msf_modoki_sum,2)
msf_modoki_atu_avg=dim_avg_n_Wrap( msf_modoki_atu,2)
msf_modoki_win_avg=dim_avg_n_Wrap( msf_modoki_win,2)


ud_modoki_sp_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modoki_month_sp),2)
ud_modoki_sum_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modoki_month_sum),2)
ud_modoki_atu_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modoki_month_atu),2)
ud_modoki_win_avg=dim_avg_n_Wrap( ud_time(:,:,ind_modoki_month_win),2)

w_modoki_sp_avg=dim_avg_n_Wrap( w_time(:,:,ind_modoki_month_sp),2)
w_modoki_sum_avg=dim_avg_n_Wrap( w_time(:,:,ind_modoki_month_sum),2)
w_modoki_atu_avg=dim_avg_n_Wrap( w_time(:,:,ind_modoki_month_atu),2)
w_modoki_win_avg=dim_avg_n_Wrap( w_time(:,:,ind_modoki_month_win),2)
   

;**********************************************************************
 ;画图
;************************************************************************
  wks=gsn_open_wks("png", ""+diri+"era_modoki四季")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="modoki_MAM"
    res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    res@cnLinePalette="BlueDarkRed18"
    res@cnLineLabelsOn=True
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
     ; res@cnLabelScalingMode="ScaleFactor"
     ; res@cnLabelScaleValueF=10^9
      
    res@lbOrientation="Vertical"
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(100,290,40)
    res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res@tmXBMinorOn=False
    res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    res@vcRefLengthF    = 0.04        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.02       ; thin out vectors
    res@vcMapDirection  = False

  
    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15

    res@vcRefAnnoOrthogonalPosF  = -1.05
    res@vcRefAnnoPerimOn=False
    res@tmYRMode ="Automatic"
  plot=new(4,graphic)

   plot(0)=gsn_csm_pres_hgt_vector(wks,msf_modoki_sp_avg,ud_modoki_sp_avg,w_modoki_sp_avg,res)
    res@gsnRightString="modoki_JJA"
   plot(1)=gsn_csm_pres_hgt_vector(wks,msf_modoki_sum_avg,ud_modoki_sum_avg,w_modoki_sum_avg,res)
    res@gsnRightString="modoki_SON"
   plot(2)=gsn_csm_pres_hgt_vector(wks,msf_modoki_atu_avg,ud_modoki_atu_avg,w_modoki_atu_avg,res)
    res@gsnRightString="modoki_DJF"
   plot(3)=gsn_csm_pres_hgt_vector(wks,msf_modoki_win_avg,ud_modoki_win_avg,w_modoki_win_avg,res)
                                                 
;-------------------------------
delete(res@gsnRightString)
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 4.
res@cnFillOn        = False
res@tmYRMode ="Automatic"

plot1=new(4,graphic)
plot1(0)= gsn_csm_pres_hgt(wks, msf_modoki_sp_avg, res)
plot1(1) = gsn_csm_pres_hgt(wks, msf_modoki_sum_avg, res)
plot1(2)= gsn_csm_pres_hgt(wks, msf_modoki_atu_avg, res)
plot1(3)= gsn_csm_pres_hgt(wks, msf_modoki_win_avg, res)
;-----------------------------------------------------
 ;去掉height字体
;----------------------------------------------------
do i=0,3
overlay(plot(i),plot1(i))
end do 
pres=True
gsn_panel(wks, plot, (/2,2/), pres)


end
