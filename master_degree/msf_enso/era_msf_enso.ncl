;**********************************************************************
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
;-----------------------------------------------------------------
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
;*****************************************************************8
begin
diri="/mnt/f/share/ersst/"
fils0=systemfunc("ls "+diri+"nino3_4.ascii")
nino3_4=asciiread(fils0, 457, "float")
nino3_4_de=dtrend_n(nino3_4(0:407), True,0)
printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)


enino=ind(nino3_4_nor.gt.1.0)
lanila=ind(nino3_4_nor.lt.-1.0)
Neutral=ind(nino3_4_nor.le.1.0.and.nino3_4_nor.ge.-1.0)

diri1="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri1+"msf_era_level_location1979-2017.nc")
f1=addfile(fils1, "r")

x=f1->locationxi_a(0:407)


x_de=dtrend_n(clm_ano(x),True, 0)
copy_VarCoords(x, x_de)

x_enino=dim_avg_n_Wrap(x_de(enino),0)
x_lanina=dim_avg_n_Wrap(x_de(lanila),0)
x_Neutral=dim_avg_n_Wrap(x_de(Neutral),0)
x_clm=dim_avg_n_Wrap(x, 0)

print("x_enino="+x_enino)
print("x_lanina="+x_lanina)
print("x_Neutral="+x_Neutral)
print("x_clm="+x_clm)







diri3="/mnt/f/share/era/"
fils3=systemfunc("ls "+diri3+"msf_time_era.nc")
f3=addfile(fils3, "r")

msf_time=f3->msf_time(:,:,0:407)
ud_time=f3->ud_time(:,:,0:407)
w_time=f3->w_time(:,:,0:407)
level=f3->level(:)
lon=f3->lon(:)


  msf_time&level@units="hpa"
  ud_time&level@units="hpa"
  w_time&level@units="hpa"
 
  msf_enino=dim_avg_n_Wrap(msf_time(:,:,enino(:)),2)
  msf_lanina=dim_avg_n_Wrap(msf_time(:,:,lanila),2)
  msf_clm=dim_avg_n_Wrap(msf_time, 2)
  msf_diff=msf_enino-msf_lanina 
  copy_VarMeta(msf_enino, msf_diff)

  msf_scale_enino=msf_enino/10^9
  msf_scale_lanina=msf_lanina/10^9
  msf_scale_clm=msf_clm/10^9
  msf_scale_diff=msf_diff/10^9

  copy_VarMeta(msf_enino, msf_scale_enino)
  copy_VarMeta(msf_lanina, msf_scale_lanina)
  copy_VarMeta(msf_clm, msf_scale_clm)
  copy_VarMeta(msf_diff, msf_scale_diff)


  ud_enino=dim_avg_n_Wrap(ud_time(:,:,enino(:)),2)
  ud_lanina=dim_avg_n_Wrap(ud_time(:,:,lanila),2)
  ud_clm=dim_avg_n_Wrap(ud_time, 2)
  ud_diff=ud_enino-ud_lanina
  printVarSummary(ud_enino)
  copy_VarMeta(ud_enino, ud_diff)
  w_enino=dim_avg_n_Wrap(w_time(:,:,enino),2)
  w_lanina=dim_avg_n_Wrap(w_time(:,:,lanila),2)
  w_clm=dim_avg_n_Wrap(w_time, 2)
  w_diff=w_enino-w_lanina
  copy_VarMeta(w_enino, w_diff)
  
  ;******************************************************************************************
  ;msf detrend场
  ;*******************************************************************************************

   msf_time_de=dtrend_n(rmMonAnnCycLLT(msf_time),True,2)
   copy_VarCoords(msf_time, msf_time_de)
   
   msf_enino_de=msf_time_de(:,:,enino)
   msf_lanina_de=msf_time_de(:,:,lanila)
   msf_neutral_de=msf_time_de(:,:,Neutral)
 
   msf_scale_enino_de=msf_enino_de/10^8
   msf_scale_lanina_de=msf_lanina_de/10^8
   msf_scale_neutral_de= msf_neutral_de/10^8

   copy_VarMeta(msf_enino_de, msf_scale_enino_de)
   copy_VarMeta(msf_lanina_de, msf_scale_lanina_de)
   copy_VarMeta(msf_neutral_de, msf_scale_neutral_de)


   enino_in_mean=dim_avg_n_Wrap(nino3_4_nor(enino),0)
   lanina_in_mean=dim_avg_n_Wrap(nino3_4_nor(lanila),0)
   neutral_in_mean=dim_avg_n_Wrap(nino3_4_nor(Neutral),0)
   dim=dimsizes(msf_scale_enino_de)
   mean1=new((/dim(0),dim(1)/),typeof(msf_scale_enino_de))
   mean1=enino_in_mean
  delete(dim)
   dim=dimsizes(msf_scale_lanina_de)
   mean2=new((/dim(0),dim(1)/),typeof(msf_scale_lanina_de))
   mean2=lanina_in_mean
   delete(dim)
   dim=dimsizes(msf_scale_neutral_de)
   mean3=new((/dim(0),dim(1)/),typeof(msf_scale_enino_de))
   mean3=neutral_in_mean   
  delete(dim)
;******************************************************************************************
;u w detrend anomoly场
;*******************************************************************************************
 ud_time_de=dtrend_n(rmMonAnnCycLLT(ud_time),True,2)
 w_time_de=dtrend_n(rmMonAnnCycLLT(w_time),True,2)

   ud_enino_de=dim_avg_n(ud_time_de(:,:,enino),2)
   ud_lanina_de=dim_avg_n(ud_time_de(:,:,lanila),2)
   ud_neutral_de=dim_avg_n(ud_time_de(:,:,Neutral),2)

   w_enino_de=dim_avg_n(w_time_de(:,:,enino),2)
   w_lanina_de=dim_avg_n( w_time_de(:,:,lanila),2)
   w_neutral_de=dim_avg_n(w_time_de(:,:,Neutral),2)
   printVarSummary(ud_enino_de)

   copy_VarCoords(ud_time_de(:,:,1), ud_enino_de)
   copy_VarCoords(ud_time_de(:,:,1), ud_lanina_de)
   copy_VarCoords(ud_time_de(:,:,1), ud_neutral_de)

   copy_VarCoords(w_time_de(:,:,1), w_enino_de)
   copy_VarCoords(w_time_de(:,:,1), w_lanina_de)
   copy_VarCoords(w_time_de(:,:,1), w_neutral_de)

   
;-------------------------------------------------------------------------
 ;用平均的enino指数做标准化
;---------------------------------------------------------------------------
;   stv1=new((dimsizes(msf_scale_enino_de)),typeof(msf_scale_enino_de))
;   stv2=new((dimsizes(msf_scale_lanina_de)),typeof(msf_scale_lanina_de))

;   do t=0,dimsizes(msf_scale_enino_de&time)-1
;      stv1(:,:,t)=(msf_scale_enino_de(:,:,t)-mean1)^2
;    end do

;   std1=sqrt(dim_sum_n(stv1, 2)/(dimsizes(msf_scale_enino_de&time)-1))
;   msf_enino_std=new((dimsizes(msf_scale_enino_de)),typeof(msf_scale_enino_de))

;   do t=0,dimsizes(msf_scale_enino_de&time)-1
;      msf_enino_std(:,:,t)=(msf_scale_enino_de(:,:,t)-mean1)/std1
;    end do
;   copy_VarCoords(msf_scale_enino_de, msf_enino_std)

; ;----------------------------------------------------------------------------
;    do t=0,dimsizes(msf_scale_lanina_de&time)-1
;      stv2(:,:,t)=(msf_scale_lanina_de(:,:,t)-mean2)^2
;    end do
;   std2=sqrt(dim_sum_n(stv2, 2)/(dimsizes(msf_scale_lanina_de&time)-1))
;   msf_lanina_std=new((dimsizes(msf_scale_lanina_de)),typeof(msf_scale_lanina_de))

;   do t=0,dimsizes(msf_scale_lanina_de&time)-1
;      msf_lanina_std(:,:,t)=(msf_scale_lanina_de(:,:,t)-mean2)/std2
;    end do
;   copy_VarCoords(msf_scale_lanina_de, msf_lanina_std)


;  msf_enino_tr=dim_avg_n_Wrap(msf_enino_std,2)
;  msf_lanina_tr=dim_avg_n_Wrap(msf_lanina_std,2)
;  msf_sum_tr=msf_enino_tr+msf_lanina_tr
;  copy_VarCoords(msf_enino_tr, msf_sum_tr)
;---------------------------------------------------------------------------------------

 msf_scale_enino_de_avg=dim_avg_n_Wrap(msf_scale_enino_de, 2)
 msf_scale_lanina_de_avg=dim_avg_n_Wrap(msf_scale_lanina_de, 2)
 msf_scale_neutral_de_avg=dim_avg_n_Wrap(msf_scale_neutral_de, 2)
printVarSummary(msf_scale_neutral_de_avg)
 ; msf_scale_diff=msf_scale_enino_de_avg+msf_scale_lanina_de_avg
 l=dimsizes(msf_scale_enino_de&level)
 m=dimsizes(msf_scale_lanina_de&lon)
 k=l*m-1
stv1=new((dimsizes(msf_scale_enino_de_avg)),typeof(msf_scale_enino_de_avg))
stv2=new((dimsizes(msf_scale_lanina_de_avg)),typeof(msf_scale_lanina_de_avg))
stv3=new((dimsizes(msf_scale_neutral_de_avg)),typeof(msf_scale_lanina_de_avg))

stv1=(msf_scale_enino_de_avg-mean1)
std1=sqrt(sum((stv1)^2)/k)
msf_enino_tr=stv1/std1 

stv2=msf_scale_lanina_de_avg - mean2
std2=sqrt(sum((stv2)^2/k))
msf_lanina_tr=stv2/std2 

stv3=msf_scale_neutral_de_avg - mean3
std3=sqrt(sum((stv3)^2/k))
msf_neutral_tr=stv3/std3 

msf_sum_tr=msf_enino_tr+msf_lanina_tr


copy_VarCoords(msf_scale_enino_de_avg, msf_enino_tr)
copy_VarCoords(msf_scale_lanina_de_avg, msf_lanina_tr) 
copy_VarCoords(msf_scale_neutral_de_avg, msf_neutral_tr)
copy_VarCoords(msf_enino_tr, msf_sum_tr)



;********************************************************************************************************
  wks1=gsn_open_wks("png", ""+diri+"msf_enso_原始场")
  
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="EI Nino"
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
    res@lbLabelBarOn=False
       
    res@vpWidthF=0.7
 
     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
      
    ; res@lbOrientation="Vertical"
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
    res@tmYRMode             = "Automatic"          ; turn off special labels on right axis
   plot=new(7,graphic)
   plot(0)=gsn_csm_pres_hgt_vector(wks1,msf_scale_enino,ud_enino,w_enino,res)
      res@gsnRightString="La Nina"
   plot(1)=gsn_csm_pres_hgt_vector(wks1,msf_scale_lanina,ud_lanina,w_lanina,res)
   res@vcRefMagnitudeF = 2.0
      res@gsnRightString="EI Nino-La Nina"
   plot(2)=gsn_csm_pres_hgt_vector(wks1,msf_scale_diff,ud_diff,w_diff,res)
   
;-------------------------------
delete(res@gsnRightString)
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 4.
res@cnFillOn        = False

plot_ol=new(3,graphic)
plot_ol(0)= gsn_csm_pres_hgt(wks1, msf_scale_enino, res)
plot_ol(1)= gsn_csm_pres_hgt(wks1, msf_scale_lanina, res)
plot_ol(2)= gsn_csm_pres_hgt(wks1, msf_scale_diff, res)

  do i=0,2
overlay(plot(i),plot_ol(i))
  end do

pres=True
pres@gsnMaximize=True
; pres@gsnPanelLabelBar=True
gsn_panel(wks1, plot(0:2), (/1,3/), pres)
delete(pres)

;----------------------------------------------------------------------
;画随时间变化
;----------------------------------------------------------------------
wks2=gsn_open_wks("png", ""+diri+"msf_enso_ano场")
    ; gsn_define_colormap(wks,"BlWhRe")
    gsn_define_colormap(wks2,"hotcolr_19lev")
    ; gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="era"
    res2@gsnRightString="1979-2017"
    ; res2@lbLabelStrings=ispan(-6,6,1)
  
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5.0
    res2@cnMaxLevelValF = 5.0
    res2@cnLevelSpacingF = 1.0
    ; res2@cnLineLabelFontHeightF=0.01
    res2@vpWidthF=0.7
    res2@cnFillOn = True 
     res2@cnMonoFillColor=False
     res2@cnFillColors=(/2,3,4,5,6,10,10,14,16,17,19,20/)
    res2@cnLineLabelsOn=False
    ; res2@cnMonoFillScale=True
    res2@cnInfoLabelOn=False

     ; res2@cnExplicitLineLabelsOn=True
     ; res2@cnLineLabelStrings=""+ispan(-2,6,1)+""
     ; res2@cnLabelScalingMode="ScaleFactor"
     ; res2@cnLabelScaleValueF=10^9
    res2@tmYRMode             = "Automatic"          ; turn off special labels on right axis
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
    res2@tmXBValues=ispan(100,290,40)
    res2@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res2@tmXBMinorOn=False

 


    plot(3)= gsn_csm_pres_hgt(wks2, msf_enino_tr, res2)
    plot(4)= gsn_csm_pres_hgt(wks2, msf_lanina_tr , res2)
    delete(res2@cnLevelSpacingF)
    delete(res2@cnLevelSelectionMode)
    delete(res2@cnMinLevelValF)
    delete(res2@cnMaxLevelValF)
    delete(res2@cnLevelSpacingF)
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5.0
    res2@cnMaxLevelValF = 5.0
    res2@cnLevelSpacingF = 1.0
    ; delete(res2@cnFillColors)
     ; res2@cnFillColors=(/3,4,5,6,10,10,14,16,17,20/)
    plot(5)= gsn_csm_pres_hgt(wks2, msf_scale_neutral_de_avg, res2)
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5.0
    res2@cnMaxLevelValF = 5.0
    res2@cnLevelSpacingF = 1.0

    plot(6)= gsn_csm_pres_hgt(wks2, msf_sum_tr, res2)

 ;--------------------------------------------------------
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLinesOn=True
  res2@cnLevels = ispan(-6,6,1)
  res2@cnLineLabelsOn=True
  res2@cnLineLabelBackgroundColor="Transparent"
  res2@cnLineThicknessF = 3.
  res2@cnFillOn        = False
  res2@gsnContourNegLineDashPattern=1
  plot_clm=new(4,graphic)

  plot_clm(0)= gsn_csm_pres_hgt(wks2, msf_scale_clm, res2) 
  plot_clm(1)= gsn_csm_pres_hgt(wks2, msf_scale_clm, res2) 
  plot_clm(2)= gsn_csm_pres_hgt(wks2, msf_scale_clm, res2) 
  plot_clm(3)= gsn_csm_pres_hgt(wks2, msf_scale_clm, res2) 

do i=3,6
  overlay(plot(i),plot_clm(i-3))
 end do
pres=True
pres@gsnMaximize=True
gsn_panel(wks2, plot(3:6), (/1,4/), pres)
delete(pres)

end
   
