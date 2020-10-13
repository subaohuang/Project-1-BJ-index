function normoalized(data,mean)
begin

msf_enino_tr=data/mean
copy_VarCoords(data, msf_enino_tr)
return msf_enino_tr
end 
begin
diri="/mnt/f/share/ersst/"
diri_pi="/mnt/f/share/era/"
f0=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
date=f0->time
lon=f0->longitude
lat=f0->latitude

time=cd_calendar(date, 1)
timeind=ind(time.le.201012.and.time.ge.195001)
sst_r=lonFlip(f0->sst(timeind,:,:))

sst=sst_r(:,{10.5:-10.5},{159.5:280.5})
; sst_area_r=sst_r(:,{30.5:-30.5},{99.5:290.5})
sst_area_r=sst_r(:,{30.5:-30.5},:)


 xo=ispan(160,280,1)
 yo=ispan(-10,10,1)

 sst_int=linint2_Wrap(sst&longitude, sst&latitude(::-1), sst(:,::-1,:), False, xo, yo, 0)

 printVarSummary(sst_int)
 ; ssta_de=dtrend_n(rmMonAnnCycLLT(sst_int(LATITUDE|:,LONGITUDE|:,time|:)),True,2)
 ssta_de=rmMonAnnCycLLT(sst_int(LATITUDE|:,LONGITUDE|:,time|:))

  nino3_4=dim_avg_n_Wrap(ssta_de({-5:5},{190:240},:), (/0,1/))
  nino4=dim_avg_n_Wrap(ssta_de({-5:5},{160:210},:), (/0,1/))
  nino1_2=dim_avg_n_Wrap(ssta_de({-10:0},{270:280},:), (/0,1/))
  nino3=dim_avg_n_Wrap(ssta_de({-5:5},{210:270},:), (/0,1/))

 opt=diri+"Hadley_nino3_4_1950_2010_nodtrend.ascii"
 asciiwrite(opt, sprintf ("%5.2f,", nino3_4))

 opt1=diri+"Hadley_nino4_1950_2010_nodtrend.ascii"
 asciiwrite(opt1, sprintf ("%5.2f,", nino4))

 opt2=diri+"Hadley_nino1_2_1950_2010_nodtrend.ascii"
 asciiwrite(opt2, sprintf ("%5.2f,", nino1_2))

 opt3=diri+"Hadley_nino3_1950_2010_nodtrend.ascii"
 asciiwrite(opt3, sprintf ("%5.2f,", nino3))

exit()

 nino3_4!0="time"
 nino3_4&time=ispan(1,732,1)
 nino3_4_nor=dim_standardize_n(nino3_4, 1, 0)


;---------------------------------------------------
 str_enino=ind(nino3_4.gt.1.0)
 str_lanila=ind(nino3_4.lt.-1.0)
 weak_enino=ind(nino3_4.gt.0.5.and.nino3_4.lt.1.0)
 weak_lanila=ind(nino3_4.gt.-1.0.and.nino3_4.lt.-0.5)
 print(str_enino)

 ssta_area=dtrend_n(rmMonAnnCycLLT(sst_area_r(latitude|:,longitude|:,time|:)),True,2)
 printVarSummary(ssta_area)
 copy_VarCoords(sst_area_r(latitude|:,longitude|:,time|:), ssta_area)

 ssta_str_enino=dim_avg_n_Wrap(ssta_area(:,:,str_enino),2)
 ssta_str_lanina=dim_avg_n_Wrap( ssta_area(:,:,str_lanila),2)
 ssta_weak_enino=dim_avg_n_Wrap(   ssta_area(:,:,weak_enino),2)
 ssta_weak_lanina=dim_avg_n_Wrap(  ssta_area(:,:,weak_lanila),2)                       
 printVarSummary(ssta_str_enino)
;-----------------------------------------------------------
;标准化
;-------------------------------------------------------------
   dim=dimsizes(ssta_str_enino)
   printVarSummary(nino3_4)
   mean1=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean1=dim_avg_n_Wrap(nino3_4(str_enino), 0)
 
   mean2=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean2=dim_avg_n_Wrap(nino3_4(str_lanila), 0)

  
   mean3=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean3=dim_avg_n_Wrap(nino3_4(weak_enino), 0) 

   mean4=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean4=dim_avg_n_Wrap(nino3_4(weak_lanila), 0)

   ssta_str_enino_nor=normoalized(ssta_str_enino,mean1)
   ssta_str_lanina_nor=normoalized(ssta_str_lanina,mean2)
   ssta_weak_enino_nor=normoalized(ssta_weak_enino,mean3)
   ssta_weak_lanina_nor=normoalized(ssta_weak_lanina,mean4)
   

   ssta_sum_str_nor=ssta_str_enino_nor-ssta_str_lanina_nor
   ssta_sum_weak_nor= ssta_weak_enino_nor -ssta_weak_lanina_nor 

  copy_VarCoords(ssta_str_enino, ssta_sum_str_nor)
  copy_VarCoords(ssta_str_lanina, ssta_sum_weak_nor)


RightString=(/"(a)","(b)","(c)","(d)","(e)","(f)"/)
MainString=(/"strong Ei nino","strong La nina","diff (a)-(b)","weak Ei nino","weak La nina","diff (e)-(d)"/)
;------------------------------------------------------------------
;画图
;-----------------------------------------------------------------
wks=gsn_open_wks("pdf", "/mnt/f/share/era/Hadley_sst_ano_enso")
; gsn_define_colormap(wks,"testcmap")
res=True
res = True
res@gsnAddCyclic=False
res@gsnDraw=False
res@gsnFrame=False
res@gsnLeftString=RightString(0)
res@tiMainString=MainString(0)
; res@cnFillColors=(/2,7,13,19,45,49,50,63,67,81,100,100,113,119,145,152,170,178,188,193,195,199/)
; res@cnFillPalette  = "temp_diff_18lev"
 ; res@cnFillPalette="BlueDarkRed18"
 ; res@cnFillPalette="BlAqGrWh2YeOrReVi22"
; res@cnFillPalette="NCV_blu_red"
; res@cnFillPalette="testcmap"
res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
res@cnFillOn = True 
res@cnLinesOn =True
res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
res@cnMinLevelValF       =   -1.0
res@cnMaxLevelValF       =   1.0
res@cnLevelSpacingF      =  0.1

res@lbOrientation = "Vertical"

res@mpLimitMode           = "LatLon"
res@mpMinLatF             = min(ssta_str_enino&latitude)
res@mpMaxLatF             = max(ssta_str_enino&latitude)
res@mpMinLonF             = min(ssta_str_enino&longitude)
res@mpMaxLonF             = max(ssta_str_enino&longitude)
res@mpCenterLonF          =avg(ssta_str_enino&longitude)

res@vpHeightF=0.4
res@vpWidthF=0.3

plot=new(6,graphic)

plot(0) = gsn_csm_contour_map(wks, ssta_str_enino_nor, res)
res@gsnLeftString=RightString(1)
res@tiMainString=MainString(1)
plot(1) = gsn_csm_contour_map(wks, ssta_str_lanina_nor, res)
res@gsnLeftString=RightString(2)
res@tiMainString=MainString(2)
plot(2) = gsn_csm_contour_map(wks, ssta_sum_str_nor, res)
res@gsnLeftString=RightString(3)
res@tiMainString=MainString(3)
plot(3) = gsn_csm_contour_map(wks, ssta_weak_enino_nor, res)
res@gsnLeftString=RightString(4)
res@tiMainString=MainString(4)
plot(4) = gsn_csm_contour_map(wks, ssta_weak_lanina_nor, res)

; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; res@cnMinLevelValF       =   -0.5
; res@cnMaxLevelValF       =   0.5
; res@cnLevelSpacingF      =  0.05

res@gsnLeftString=RightString(5)
res@tiMainString=MainString(5)
plot(5) = gsn_csm_contour_map(wks, ssta_sum_weak_nor, res)
printVarSummary(plot(0))
j=new(1,graphic)
j=plot(0)
 getvalues j@contour                                                     
  "cnFillColors" : am_ids                                     
  end getvalues           
 k=am_ids
 opt=diri+"rgb.ascii"
 asciiwrite(opt, sprintf ("%10.0f,", k))
pres=True
pres@gsnPanelDebug=True

gsn_panel(wks, plot, (/2,3/), pres)


end




