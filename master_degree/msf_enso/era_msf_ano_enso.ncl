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
;-------------------------------------
function normoalized(data,mean)
begin

msf_enino_tr=data/mean
copy_VarCoords(data, msf_enino_tr)
msf_enino_tr!0="lev"
msf_enino_tr&lev@units = "hpa"
return msf_enino_tr
end 
;---------------------------------------

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
;-------------------------------------------------

begin
diri="/mnt/f/share/era/data/"
diri_pi="/mnt/f/share/era/"
fils0=systemfunc("ls "+diri+"era_nino3_4.ascii")
nino3_4=asciiread(fils0, 456, "float")
nino3_4_de=dtrend_n(nino3_4(0:407), True,0)
printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)

 enino=ind(nino3_4_de.gt.1.0)
 lanila=ind(nino3_4_de.lt.-1.0)
 ; Neutral=ind(nino3_4_de.lt.1.0.and.nino3_4_de.gt.-1.0)

 ; str_enino=ind(nino3_4_nor.gt.1.0)
 ; str_lanila=ind(nino3_4_nor.lt.-1.0)
 ; weak_enino=ind(nino3_4_nor.gt.0.5.and.nino3_4_nor.lt.1.0)
 ; weak_lanila=ind(nino3_4_nor.gt.-1.0.and.nino3_4_nor.lt.-0.5)

 str_enino=ind(nino3_4_de.gt.1.0)
 str_lanila=ind(nino3_4_de.lt.-1.0)
 weak_enino=ind(nino3_4_de.gt.0.5.and.nino3_4_de.lt.1.0)
 weak_lanila=ind(nino3_4_de.gt.-1.0.and.nino3_4_de.lt.-0.5)
 Neutral=ind(nino3_4_de.lt.0.5.and.nino3_4_de.gt.-0.5)
  

diri1="/mnt/f/share/era/data/"

; fils1=systemfunc("ls "+diri1+"msf_era_level_location1979-2017.nc")
fils1=systemfunc("ls "+diri1+"msf_era_level_location1979-2017_reconsitution.nc")
f1=addfile(fils1, "r")
x=f1->locationxi_a(0:407)
x_de=dtrend_n(clm_ano(x), True, 0)
print(x_de)







diri3="/mnt/f/share/era/data/"
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

   msf_str_enino_de=msf_time_de(:,:,str_enino)
   msf_str_lanina_de=msf_time_de(:,:,str_lanila)
   msf_weak_enino_de=msf_time_de(:,:,weak_enino)
   msf_weak_lanina_de=msf_time_de(:,:,weak_lanila)
 
   msf_str_scale_enino_de=msf_str_enino_de/10^9
   msf_str_scale_lanina_de=msf_str_lanina_de/10^9
   msf_weak_scale_enino_de=msf_weak_enino_de/10^9
   msf_weak_scale_lanina_de=msf_weak_lanina_de/10^9
   
   msf_scale_time_de= msf_time_de/10^9
   
   copy_VarMeta(msf_str_enino_de, msf_str_scale_enino_de)
   copy_VarMeta(msf_str_lanina_de, msf_str_scale_lanina_de)
   copy_VarMeta(msf_weak_enino_de, msf_weak_scale_enino_de)
   copy_VarMeta(msf_weak_lanina_de, msf_weak_scale_lanina_de)

   copy_VarMeta(msf_time_de, msf_scale_time_de)


   str_enino_in_mean=dim_avg_n_Wrap(nino3_4_nor(str_enino),0)
   str_lanina_in_mean=dim_avg_n_Wrap(nino3_4_nor(str_lanila),0)
   weak_enino_in_mean=dim_avg_n_Wrap(nino3_4_nor(weak_enino),0)
   weak_lanina_in_mean=dim_avg_n_Wrap(nino3_4_nor(weak_lanila),0)

   clm_in_mean=dim_avg_n_Wrap(nino3_4_nor, 0)
  
   dim=dimsizes(msf_str_scale_enino_de)
   mean1=new((/dim(0),dim(1)/),typeof(msf_str_scale_enino_de))
   mean1=str_enino_in_mean
  delete(dim)
   dim=dimsizes(msf_str_scale_enino_de)
   mean2=new((/dim(0),dim(1)/),typeof(msf_str_scale_enino_de))
   mean2=str_lanina_in_mean
   delete(dim)
  dim=dimsizes(msf_str_scale_enino_de)
   mean3=new((/dim(0),dim(1)/),typeof(msf_str_scale_enino_de))
   mean3=weak_enino_in_mean
  delete(dim)
 dim=dimsizes(msf_str_scale_enino_de)
   mean4=new((/dim(0),dim(1)/),typeof(msf_str_scale_enino_de))
   mean4=weak_lanina_in_mean
  delete(dim)
;-----------------------------------------------------------------------

x_str_enino=dim_avg_n_Wrap(x_de(str_enino),0)
x_str_lanina=dim_avg_n_Wrap(x_de(str_lanila),0)
x_Neutral=dim_avg_n_Wrap(x_de(Neutral),0)
x_weak_enino=dim_avg_n_Wrap(x_de(weak_enino),0)
x_weak_lanina=dim_avg_n_Wrap(x_de(weak_lanila),0)
x_clm=dim_avg_n_Wrap(x_de,0)


print("x_str_enino="+x_str_enino/str_enino_in_mean)
print("x_str_lanina="+x_str_lanina/str_lanina_in_mean)
print("x_weak_enino="+x_weak_enino/weak_enino_in_mean)
print("x_weak_lanina="+x_weak_lanina/weak_lanina_in_mean)
print("x_Neutral="+x_Neutral)
print("x_clm="+x_clm)

; print("x_str_enino="+x_str_enino)
; print("x_str_lanina="+x_str_lanina)
; print("x_weak_enino="+x_weak_enino)
; print("x_weak_lanina="+x_weak_lanina)
; print("x_Neutral="+x_Neutral)
; print("x_clm="+x_clm)


 





;------------------------------------------------------------------------
;ud detrend anomoly场
;--------------------------------------------------------------------------
   ud_time_de=dtrend_n(rmMonAnnCycLLT(ud_time),True,2)
   w_time_de=dtrend_n(rmMonAnnCycLLT(w_time),True,2)
   copy_VarCoords(ud_time, ud_time_de)
   copy_VarCoords(w_time, w_time_de)

   ud_enino_de=dim_avg_n_Wrap(ud_time_de(:,:,enino),2)
   ud_lanina_de=dim_avg_n_Wrap( ud_time_de(:,:,lanila),2)
   ud_neutral_de=dim_avg_n_Wrap(ud_time_de(:,:,Neutral),2)
   ud_clm_de=dim_avg_n_Wrap(ud_time, 2)


   ud_str_enino_de=dim_avg_n_Wrap(ud_time_de(:,:,str_enino),2)
   ud_str_lanina_de=dim_avg_n_Wrap(ud_time_de(:,:,str_lanila),2)
   
   ud_weak_enino_de=dim_avg_n_Wrap(ud_time_de(:,:,weak_enino),2)
   ud_weak_lanina_de=dim_avg_n_Wrap(ud_time_de(:,:,weak_lanila),2)
   
   printVarSummary(ud_time)


   ud_neutral_de=dim_avg_n_Wrap(ud_time_de(:,:,Neutral),2)
   ud_clm_de=dim_avg_n_Wrap(ud_time_de, 2)
    
   ud_diff=ud_enino_de+ud_lanina_de
   copy_VarCoords(ud_enino_de, ud_diff)

   ud_str_diff=ud_str_enino_de- ud_str_lanina_de
   copy_VarCoords(ud_str_enino_de, ud_str_diff)
   
     
   ud_weak_diff=ud_weak_enino_de- ud_weak_lanina_de
   copy_VarCoords(ud_weak_enino_de, ud_weak_diff)





   w_enino_de=dim_avg_n_Wrap(w_time_de(:,:,enino),2)
   w_lanina_de=dim_avg_n_Wrap( w_time_de(:,:,lanila),2)
   w_neutral_de=dim_avg_n_Wrap(w_time_de(:,:,Neutral),2)
   w_clm_de=dim_avg_n_Wrap(w_time, 2)

   w_str_enino_de=dim_avg_n_Wrap(w_time_de(:,:,str_enino),2)
   w_str_lanina_de=dim_avg_n_Wrap(w_time_de(:,:,str_lanila),2)

   w_weak_enino_de=dim_avg_n_Wrap(w_time_de(:,:,weak_enino),2)
   w_weak_lanina_de=dim_avg_n_Wrap(w_time_de(:,:,weak_lanila),2)
    
   w_diff=w_enino_de+w_lanina_de
   copy_VarCoords(w_enino_de, w_diff)



   w_diff=w_enino_de+w_lanina_de
   copy_VarCoords(w_enino_de, w_diff)
   
  

   
  ud_str_enino_tr=normoalized(ud_str_enino_de,mean1)
  ud_str_lanina_tr=normoalized(ud_str_lanina_de,mean2)
  ud_weak_enino_tr=normoalized(ud_weak_enino_de,mean3)
  ud_weak_lanina_tr=normoalized(ud_weak_lanina_de,mean4)
  

  w_str_enino_tr=normoalized(w_str_enino_de,mean1)
  w_str_lanina_tr=normoalized(w_str_lanina_de,mean2)
  w_weak_enino_tr=normoalized(w_weak_enino_de,mean3)
  w_weak_lanina_tr=normoalized(w_weak_lanina_de,mean4)
  

   w_str_diff=w_str_enino_tr- w_str_lanina_tr
   copy_VarCoords(w_str_enino_tr, w_str_diff)
   
     
   w_weak_diff=w_weak_enino_tr- w_weak_lanina_tr

   copy_VarCoords(w_weak_enino_tr, w_weak_diff)




;---------------------------------------------------------------------------------------
;用平均的mean index相除做norm化。
;--------------------------------------------------------------------------------------

  msf_str_enino_tr=normoalized(dim_avg_n_Wrap(msf_str_scale_enino_de,2),mean1)
  msf_str_lanina_tr=normoalized(dim_avg_n_Wrap(msf_str_scale_lanina_de,2),mean2)
  msf_weak_enino_tr=normoalized(dim_avg_n_Wrap( msf_weak_scale_enino_de,2),mean3)
  msf_weak_lanina_tr=normoalized(dim_avg_n_Wrap( msf_weak_scale_lanina_de,2),mean4)
  
 
  msf_str_enino_clm=dim_avg_n_Wrap(msf_str_scale_enino_de,2)
  msf_weak_enino_clm=dim_avg_n_Wrap(msf_weak_scale_enino_de,2)
  msf_str_enino_clm=where(msf_str_enino_clm.lt.0, msf_str_enino_clm, msf_str_enino_clm@_FillValue)
  msf_weak_enino_clm=where(msf_weak_enino_clm.lt.0, msf_weak_enino_clm, msf_weak_enino_clm@_FillValue)
 
  msf_str_lanina_clm=dim_avg_n_Wrap(msf_str_scale_lanina_de,2)
  msf_weak_lanina_clm=dim_avg_n_Wrap(msf_weak_scale_lanina_de,2)
  msf_str_lanina_clm=where(msf_str_lanina_clm.gt.0, msf_str_lanina_clm, msf_str_lanina_clm@_FillValue)
  msf_weak_lanina_clm=where(msf_weak_lanina_clm.gt.0, msf_weak_lanina_clm, msf_weak_lanina_clm@_FillValue)



  str_enino_clm=avg(msf_str_enino_clm)
  print("str_enino_clm="+str_enino_clm)

  weak_enino_clm=avg(msf_weak_enino_clm)
  print("weak_enino_clm="+weak_enino_clm)

  str_lanina_clm=avg(msf_str_lanina_clm)
  print("str_lanina_clm="+str_lanina_clm)

  weak_lanina_clm=avg(msf_weak_lanina_clm)
  print("weak_lanina_clm="+weak_lanina_clm)





  msf_str_sum=msf_str_enino_tr - msf_str_lanina_tr
  msf_weak_sum=  msf_weak_lanina_tr-msf_weak_enino_tr

  printVarSummary(msf_str_sum)
  copy_VarCoords(msf_str_enino_tr, msf_str_sum)
  copy_VarCoords(msf_weak_enino_tr, msf_weak_sum)
;********************************************************************************************************
  wks=gsn_open_wks("png", ""+diri_pi+"era_msf_enso_去偏态")
  plot=new(6,graphic)
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
    res2@gsnRightString="strong Ei Nino"
    ; res2@lbLabelStrings=ispan(-6,6,1)
    
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -2.0
    res2@cnMaxLevelValF = 2.0
    res2@cnLevelSpacingF = 0.5

    res2@lbLabelBarOn=False
    res2@cnLineLabelFontHeightF=0.01
    res2@vpWidthF=0.7
    res2@cnFillOn = True 
    res2@cnMonoFillColor=False
    res2@cnFillColors=(/2,4,5,6,10,10,14,16,17,20/)
    res2@cnLineLabelsOn=False
    ; res2@cnMonoFillScale=True
    res2@cnInfoLabelOn=False
    res2@lbTitleString="10~S~9~N~kgs-1"
    res2@lbTitlePosition="Bottom"
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
    plot(0)= gsn_csm_pres_hgt_vector(wks, msf_str_enino_tr,ud_str_enino_tr,w_str_enino_tr ,res2)
    res2@gsnRightString="strong La Nina"
    plot(1)= gsn_csm_pres_hgt_vector(wks, msf_str_lanina_tr ,ud_str_lanina_tr,w_str_lanina_tr, res2)
    
   
    
    res2@gsnRightString="strong Ei Nino-La Nina"

    res2@cnMinLevelValF = -1.0
    res2@cnMaxLevelValF = 1.0
    res2@cnLevelSpacingF = 0.25
    delete(res2@vcRefMagnitudeF)
  
     ; res2@vpYF =0.785195
    plot(2)= gsn_csm_pres_hgt_vector(wks, msf_str_sum ,ud_str_diff,w_str_diff, res2)


  
    res2@cnMinLevelValF = -2.0
    res2@cnMaxLevelValF = 2.0
    res2@cnLevelSpacingF = 0.5
    res2@gsnRightString="weak Ei Nino "
    res2@vcRefMagnitudeF = 1.5            ; define vector ref mag
    ; delete(res2@vpWidthF)
    ; delete(res2@vpHeightF)
    ; delete(res2@vpYF)
   plot(3)= gsn_csm_pres_hgt_vector(wks, msf_weak_enino_tr ,ud_weak_enino_tr,w_weak_enino_tr, res2)
    
    res2@gsnRightString="weak La Nina"
   plot(4)= gsn_csm_pres_hgt_vector(wks, msf_weak_lanina_tr ,ud_weak_lanina_tr,w_weak_lanina_tr, res2)
    ; res2@vcRefMagnitudeF = 1.0
    delete(res2@vcRefMagnitudeF)
    res2@gsnRightString="weak Ei Nino-weak La Nina"
    
    res2@cnMinLevelValF = -1.0
    res2@cnMaxLevelValF = 1.0
    res2@cnLevelSpacingF = 0.25
   plot(5)= gsn_csm_pres_hgt_vector(wks, msf_weak_sum , ud_weak_diff,w_weak_diff, res2)

  ;---Resources for the polyline  添加下垫面的陆地
  ployres                  = True          ; polyline mods desired
  ployres@gsLineThicknessF = 15.0            ; line thickness   
  ployres@gsLineColor      = "black"         ; color of lines
  ployres@gsLineDashPattern=15
  ;---Attach the polyline
  ploylevel=1000
  ploylon1=fspan(100, 149, 25)
  ploylon2=fspan(280, 290, 11)
  ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
  ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
  
  dum1=new(6,graphic)
  dum2=new(6,graphic) 

   do i=0,5
 dum1(i)= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
  end do
ployres@gsLineDashPattern=0
do i=0,5
 dum2(i)= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel1,  ployres)
  end do

 ;--------------------------------------------------------
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLinesOn=True
  res2@cnLevels = ispan(-6,6,1)
  res2@cnLineLabelsOn=True
  res2@cnLineLabelBackgroundColor="Transparent"
  res2@cnLineThicknessF = 3.
  res2@cnFillOn        = False
  res2@gsnContourNegLineDashPattern=1
  delete(res2@gsnRightString)
  plot_clm=new(6,graphic)
do i=0,5
  plot_clm(i)= gsn_csm_pres_hgt(wks, msf_scale_clm, res2) 
end do
  
do i=0,5
  overlay(plot(i),plot_clm(i))
 end do

; pres=True
; pres@gsnMaximize=True
; ; pres@gsnPanelLabelBar=True
; gsn_panel(wks, plot, (/2,3/), pres)

 pres1                  = True                ; mods desired
 pres1@gsnFrame         = False               ; save panel until both ready
 pres1@gsnPanelRight    = 0.60               ; draw from left edge to center
 pres1@gsnPanelLabelBar = True                ; commom label bar
 pres1@gsnPanelDebug=True
 pres1@lbTitleString    = "10~S~9~N~kg~N~s~S~-1~N~k~S~-1"                ; title string
 pres1@lbTitlePosition  = "Right"              ; title position
 pres1@lbTitleFontHeightF= .015                ; make title smaller
 pres1@lbTitleDirection = "Across"             ; title direction
 pres1@gsnPanelFigureStrings=(/"a","b","d","e"/)
 pres1@gsnPanelFigureStringsJust="TopLeft"


 pres2                  = True                ; mods desired
 pres2@gsnFrame         = False               ; save panel until both ready
 pres2@gsnPanelLeft     = 0.50      ; draw from center to right edge 
 pres2@gsnPanelLabelBar = True          ; commom label bar
  pres2@gsnPanelTop=0.78
  pres2@gsnPanelBottom=0.22
  pres2@pmLabelBarWidthF=0.23
  ; pres2@pmLabelBarHeightF=0.1

  pres2@lbTitleOn        =  True                ; turn on title
 pres2@lbTitleString    = "10~S~9~N~kg~N~s~S~-1~N~k~S~-1"                ; title string
 pres2@lbTitlePosition  = "Right"              ; title position
 pres2@lbTitleFontHeightF= .015                ; make title smaller
 pres2@lbTitleDirection = "Across"             ; title direction
 pres2@gsnPanelFigureStrings=(/"c","f"/)
  pres2@gsnPanelFigureStringsJust="TopLeft"                                     
  
 gsn_panel(wks,(/plot(0),plot(1),plot(3),plot(4)/),(/2,2/),pres1)          ; create first panel
 gsn_panel(wks,(/plot(2),plot(5)/),(/2,1/),pres2)    

 frame(wks)                                   ; draw both on same page





end
   
