function normoalized(data,mean)
begin

msf_enino_tr=data/mean
copy_VarCoords(data, msf_enino_tr)
return msf_enino_tr
end 

begin
diri="/mnt/f/share/ersst/"
fils0=systemfunc("ls "+diri+"nino3_4_1950_2010.ascii")
nino3_4_de=asciiread(fils0, 732, "float")
; nino3_4_de=dtrend_n(nino3_4(0:731), True,0)

nino3_4_de!0="time"
nino3_4_de&time=ispan(1,732,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)
;---------------------------------------------------
 str_enino=ind(nino3_4_de.gt.1.0)
 str_lanila=ind(nino3_4_de.lt.-1.0)
 weak_enino=ind(nino3_4_de.gt.0.5.and.nino3_4_de.lt.1.0)
 weak_lanila=ind(nino3_4_de.gt.-1.0.and.nino3_4_de.lt.-0.5)


m=""+ispan(1950,2010,1)+""
;  fils1=systemfunc("ls "+diri+"MMEAN0089*.nc")
;**************************************************
 ;批量读取nc文件，by lin
;*************************************************
yr=dimsizes(m)
fils=new((/yr,12/),string)
do i=0,yr-1
fils1=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
fils(i,:)= fils1
end do
fils_new=onedtond(fils, yr*12)
; fileh=new(yr*12,string)
; fileh(456)="/mnt/f/share/ersst/ersst.201701.nc" 
; fileh(0:455)=fils_new
f1=addfiles(fils_new,"join")
ssta=f1[:]->ssta
printVarSummary(ssta)


 f2=addfile(fils_new(456),"r")
 lon=f2->lon(:)
 lat=f2->lat(:)
 latind1=ind(lat.eq.-30)
 latind2=ind(lat.eq.30)
 ; lonind1=ind(lon.eq.100) ;;画全球的，把这里去掉了
 ; lonind2=ind(lon.eq.290) ;;

 ssta_area=dtrend_n(ssta(:,0,latind1:latind2,:),True,0)
 copy_VarCoords(ssta(:,0,latind1:latind2,:), ssta_area)
 printVarSummary(ssta_area)

 ssta_str_enino=dim_avg_n_Wrap(ssta_area(str_enino,:,:),0)
 ssta_str_lanina=dim_avg_n_Wrap( ssta_area(str_lanila,:,:),0)
 ssta_weak_enino=dim_avg_n_Wrap(   ssta_area(weak_enino,:,:),0)
 ssta_weak_lanina=dim_avg_n_Wrap(  ssta_area(weak_lanila,:,:),0)

 ssta_sum_str=ssta_str_enino+ssta_str_lanina
 ssta_sum_weak=ssta_weak_enino+ssta_weak_lanina

 copy_VarCoords(ssta_str_enino, ssta_sum_str)
 copy_VarCoords(ssta_weak_enino, ssta_sum_weak)


 printVarSummary(ssta_str_enino)

;-----------------------------------------------------------
;标准化
;-------------------------------------------------------------
   dim=dimsizes(ssta_str_enino)
   printVarSummary(nino3_4_de)
   mean1=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean1=dim_avg_n_Wrap(nino3_4_de(str_enino), 0)
 
   mean2=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean2=dim_avg_n_Wrap(nino3_4_de(str_lanila), 0)

  
   mean3=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean3=dim_avg_n_Wrap(nino3_4_de(weak_enino), 0) 

   mean4=new((/dim(0),dim(1)/),typeof(ssta_str_enino))
   mean4=dim_avg_n_Wrap(nino3_4_de(weak_lanila), 0)

   ssta_str_enino_nor=normoalized(ssta_str_enino,mean1)
   ssta_str_lanina_nor=normoalized(ssta_str_lanina,mean2)
   ssta_weak_enino_nor=normoalized(ssta_weak_enino,mean3)
   ssta_weak_lanina_nor=normoalized(ssta_weak_lanina,mean4)
   


   ssta_sum_str_nor=ssta_str_enino_nor-ssta_str_lanina_nor
   ssta_sum_weak_nor=  ssta_weak_enino_nor - ssta_weak_lanina_nor

  copy_VarCoords(ssta_str_enino, ssta_sum_str_nor)
  copy_VarCoords(ssta_str_lanina, ssta_sum_weak_nor)


RightString=(/"(a)","(b)","(c)","(d)","(e)","(f)"/)
MainString=(/"strong Ei nino","strong La nina","diff (a)-(b)","weak Ei nino","weak La nina","diff (e)-(d)"/)
;------------------------------------------------------------------
;画图
;-----------------------------------------------------------------
 wks=gsn_open_wks("pdf", "/mnt/f/share/era/sst_ano_enso")
res=True
res = True
res@gsnAddCyclic=False
res@gsnDraw=False
res@gsnFrame=False
res@gsnLeftString=RightString(0)
res@tiMainString=MainString(0)
; res@cnFillColors=(/2,7,13,19,45,49,50,63,67,81,100,100,113,119,145,152,170,178,188,193,195,200/)
; res@cnFillPalette  = "temp_diff_18lev"
 ; res@cnFillPalette="BlueDarkRed18"
 ; res@cnFillPalette="BlAqGrWh2YeOrReVi22"
res@cnFillPalette="NCV_blu_red"
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
res@mpMinLatF             = min(ssta_str_enino&lat)
res@mpMaxLatF             = max(ssta_str_enino&lat)
res@mpMinLonF             = min(ssta_str_enino&lon)
res@mpMaxLonF             = max(ssta_str_enino&lon)
res@mpCenterLonF          =avg(ssta_str_enino&lon)

; res@vpHeightF=0.4
; res@vpWidthF=0.3

plot=new(6,graphic)

plot(0) = gsn_csm_contour_map(wks, ssta_str_enino, res)
res@gsnLeftString=RightString(1)
res@tiMainString=MainString(1)
plot(1) = gsn_csm_contour_map(wks, ssta_str_lanina, res)
res@gsnLeftString=RightString(2)
res@tiMainString=MainString(2)
plot(2) = gsn_csm_contour_map(wks, ssta_sum_str, res)
res@gsnLeftString=RightString(3)
res@tiMainString=MainString(3)
plot(3) = gsn_csm_contour_map(wks, ssta_weak_enino, res)
res@gsnLeftString=RightString(4)
res@tiMainString=MainString(4)
plot(4) = gsn_csm_contour_map(wks, ssta_weak_lanina, res)
res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
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
pres@gsnMaximize=True
gsn_panel(wks, plot, (/2,3/), pres)


end


