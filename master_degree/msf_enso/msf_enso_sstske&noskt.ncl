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

;*****************************************************************
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
fils1=systemfunc("ls "+diri1+"msf_era_gravity_center_location1979-2017.nc")
f1=addfile(fils1, "r")
x=f1->x(0:407)
printVarSummary(x)


x_enino=dim_avg_n_Wrap(x(enino),0)
x_lanina=dim_avg_n_Wrap(x(lanila),0)
x_Neutral=dim_avg_n_Wrap(x(Neutral),0)
print("x_enino="+x_enino)
print("x_lanina="+x_lanina)
print("x_Neutral="+x_Neutral)

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
   msf_scale_time_de= msf_time_de/10^8
   

   copy_VarMeta(msf_enino_de, msf_scale_enino_de)
   copy_VarMeta(msf_lanina_de, msf_scale_lanina_de)
   copy_VarMeta(msf_neutral_de, msf_scale_neutral_de)
   copy_VarMeta(msf_time_de, msf_scale_time_de)


   enino_in_mean=dim_avg_n_Wrap(nino3_4_nor(enino),0)
   lanina_in_mean=dim_avg_n_Wrap(nino3_4_nor(lanila),0)
   neutral_in_mean=dim_avg_n_Wrap(nino3_4_nor(Neutral),0)
   clm_in_mean=dim_avg_n_Wrap(nino3_4_nor, 0)
  
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
 dim=dimsizes(msf_scale_time_de)
   mean4=new((/dim(0),dim(1)/),typeof(msf_scale_enino_de))
   mean4=clm_in_mean
  delete(dim)

;------------------------------------------------------------------------
;ud detrend anomoly场
;--------------------------------------------------------------------------
   ud_time_de=dtrend_n(rmMonAnnCycLLT(ud_time),True,2)
   w_time_de=dtrend_n(rmMonAnnCycLLT(w_time),True,2)
   

   ud_enino_de=dim_avg_n(ud_time_de(:,:,enino),2)
   ud_lanina_de=dim_avg_n(ud_time_de(:,:,lanila),2)
   ud_neutral_de=dim_avg_n(ud_time_de(:,:,Neutral),2)
   ud_clm_de=dim_avg_n(ud_time_de, 2)

   ud_diff=ud_enino_de+ud_lanina_de
   copy_VarCoords(ud_enino_de, ud_diff)
   

   w_enino_de=dim_avg_n(w_time_de(:,:,enino),2)
   w_lanina_de=dim_avg_n( w_time_de(:,:,lanila),2)
   w_neutral_de=dim_avg_n(w_time_de(:,:,Neutral),2)
   w_clm_de=dim_avg_n(w_time, 2)

   w_diff=w_enino_de+w_lanina_de
   copy_VarCoords(w_enino_de, w_diff)

   copy_VarCoords(ud_time(:,:,1), ud_enino_de)
   copy_VarCoords(ud_time(:,:,1), ud_lanina_de)
   copy_VarCoords(ud_time(:,:,1), ud_neutral_de)
   copy_VarCoords(ud_time(:,:,1), ud_clm_de)
   
   copy_VarCoords(w_time(:,:,1), w_enino_de)
   copy_VarCoords(w_time(:,:,1), w_lanina_de)
   copy_VarCoords(w_time(:,:,1), w_neutral_de)
   copy_VarCoords(w_time(:,:,1), w_clm_de)
 

;---------------------------------------------------------------------------------------
;用平均的mean index做标准化
;--------------------------------------------------------------------------------------

 msf_scale_enino_de_avg=dim_avg_n_Wrap(msf_scale_enino_de, 2)
 msf_scale_lanina_de_avg=dim_avg_n_Wrap(msf_scale_lanina_de, 2)
 msf_scale_neutral_de_avg=dim_avg_n_Wrap(msf_scale_neutral_de, 2)
 msf_scale_sum_de=msf_scale_enino_de_avg+msf_scale_lanina_de_avg
 copy_VarCoords(msf_scale_enino_de, msf_scale_sum_de)
 msf_scale_time_de_avg=dim_avg_n_Wrap(msf_scale_time_de, 2)


 ; msf_scale_diff=msf_scale_enino_de_avg+msf_scale_lanina_de_avg
 l=dimsizes(msf_scale_enino_de&level)
 m=dimsizes(msf_scale_lanina_de&lon)
 k=l*m-1
stv1=new((dimsizes(msf_scale_enino_de_avg)),typeof(msf_scale_enino_de_avg))
stv2=new((dimsizes(msf_scale_lanina_de_avg)),typeof(msf_scale_lanina_de_avg))
stv3=new((dimsizes(msf_scale_neutral_de_avg)),typeof(msf_scale_lanina_de_avg))
stv4=new((dimsizes(msf_scale_time_de_avg)),typeof(msf_scale_time_de_avg))


stv1=(msf_scale_enino_de_avg-mean1)
std1=sqrt(sum((stv1)^2)/k)
msf_enino_tr=stv1/std1 

stv2=msf_scale_lanina_de_avg - mean2
std2=sqrt(sum((stv2)^2/k))
msf_lanina_tr=stv2/std2 

stv3=msf_scale_neutral_de_avg - mean3
std3=sqrt(sum((stv3)^2/k))
msf_neutral_tr=stv3/std3 

stv4=msf_scale_time_de_avg - mean4
std4=sqrt(sum((stv4)^2/k))
msf_time_tr=stv4/std4 


msf_sum_tr=msf_enino_tr+msf_lanina_tr

copy_VarCoords(msf_scale_enino_de_avg, msf_enino_tr)
copy_VarCoords(msf_scale_lanina_de_avg, msf_lanina_tr) 
copy_VarCoords(msf_scale_neutral_de_avg, msf_neutral_tr)
copy_VarCoords(msf_scale_time_de_avg, msf_time_tr)
copy_VarCoords(msf_enino_tr, msf_sum_tr)


;********************************************************************************************************
  wks=gsn_open_wks("pdf", ""+diri+"msf_enso_去偏态和不去偏态比较")
  plot=new(8,graphic)
;----------------------------------------------------------------------
;画随时间变化
;----------------------------------------------------------------------
    ; gsn_define_colormap(wks,"BlWhRe")
    gsn_define_colormap(wks,"hotcolr_19lev")
    ; gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="era"
    res2@gsnRightString="1979-2017"
    ; res2@lbLabelStrings=ispan(-6,6,1)
  
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5
    res2@cnMaxLevelValF = 5
    res2@cnLevelSpacingF = 1
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

    res2@vcRefMagnitudeF = 1.5              ; define vector ref mag
    res2@vcRefLengthF    = 0.04        ; define length of vec ref
    res2@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res2@vcMinDistanceF  = 0.02       ; thin out vectors
    res2@vcMapDirection  = False

    ; res@vcFillArrowHeadXF=0.1
    res2@vcFillArrowMinFracWidthF=1.0
    res2@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res2@vcFillArrowWidthF=0.05
    res2@vcFillArrowHeadYF=0.15

    res2@vcRefAnnoOrthogonalPosF  = -1.05
    res2@vcRefAnnoPerimOn=False
  
   ;画的是扣掉SST偏态的场。
    plot(0)= gsn_csm_pres_hgt_vector(wks, msf_enino_tr,ud_enino_de,w_enino_de, res2)
    plot(1)= gsn_csm_pres_hgt_vector(wks, msf_lanina_tr ,ud_lanina_de,w_lanina_de, res2)

    delete(res2@cnLevelSelectionMode )
    delete(res2@cnMinLevelValF )
    delete(res2@cnMaxLevelValF )
    delete(res2@cnLevelSpacingF )
    delete(res2@cnFillColors)
 res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -32
    res2@cnMaxLevelValF = 32
    res2@cnLevelSpacingF = 8
    res2@cnFillColors=(/3,4,5,6,10,10,14,16,17,19/)
    plot(4)=   gsn_csm_pres_hgt_vector(wks, msf_scale_enino_de_avg,ud_enino_de,w_enino_de, res2)
    plot(5)=   gsn_csm_pres_hgt_vector(wks, msf_scale_lanina_de_avg,ud_lanina_de,w_lanina_de, res2)
    

    delete(res2@cnLevelSpacingF)
    ; res2@cnLevelSelectionMode ="ManualLevels"
    ; res2@cnMinLevelValF = -1
    ; res2@cnMaxLevelValF = 1
    ; res2@cnLevelSpacingF = 1
     delete(res2@cnLevelSelectionMode )
    delete(res2@cnMinLevelValF )
    delete(res2@cnMaxLevelValF )
    delete(res2@cnLevelSpacingF )
    delete(res2@cnFillColors)
    ; delete(res2@cnFillColors)
     ; res2@cnFillColors=(/3,4,5,6,10,10,14,16,17,20/)
     res2@vcRefMagnitudeF = 0.15
      res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5
    res2@cnMaxLevelValF = 5
    res2@cnLevelSpacingF = 1
    res2@cnFillColors=(/2,3,4,5,6,10,10,14,16,17,19,20/)
    plot(2)= gsn_csm_pres_hgt_vector(wks, msf_neutral_tr,ud_neutral_de,w_neutral_de,res2)
    plot(6)= gsn_csm_pres_hgt_vector(wks, msf_scale_neutral_de_avg,ud_neutral_de,w_neutral_de,res2)
    res2@vcRefMagnitudeF = 2.0
  
    delete(res2@cnLevelSelectionMode )
    delete(res2@cnMinLevelValF )
    delete(res2@cnMaxLevelValF )
    delete(res2@cnLevelSpacingF )
 
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -5
    res2@cnMaxLevelValF = 5
    res2@cnLevelSpacingF = 1
   delete(res2@cnFillColors)
    plot(3)= gsn_csm_pres_hgt_vector(wks, msf_sum_tr,ud_diff,w_diff,res2)
    res2@cnMinLevelValF = -32
    res2@cnMaxLevelValF = 32
    res2@cnLevelSpacingF = 8
    res2@cnFillColors=(/3,4,5,6,10,10,14,16,17,19/)
    plot(7)= gsn_csm_pres_hgt_vector(wks, msf_scale_sum_de,ud_diff,w_diff,res2)

 ;--------------------------------------------------------
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLinesOn=True
  res2@cnLevels = ispan(-6,6,1)
  res2@cnLineLabelsOn=True
  res2@cnLineLabelBackgroundColor="Transparent"
  res2@cnLineThicknessF = 3.
  res2@cnFillOn        = False
  res2@gsnContourNegLineDashPattern=1
  plot_clm=new(8,graphic)
do i=0,7
  plot_clm(i)= gsn_csm_pres_hgt(wks, msf_scale_clm, res2) 
end do
  

  

do i=0,7
  overlay(plot(i),plot_clm(i))
 end do

pres=True
pres@gsnMaximize=True
; pres@gsnPanelLabelBar=True
gsn_panel(wks, plot, (/2,4/), pres)

end
   
