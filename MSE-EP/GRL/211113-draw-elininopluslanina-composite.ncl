
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
     ; var_pv_ano_east    =  dim_variance_n(pv_ano, 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)
    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

function mse_diff(data,year1,year2,year3,year4)
begin
  ; rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  ; copy_VarCoords(data(1,:,:), rc1)
  ; rc1 = rc1*(year2-year1+1)
  rc2 = dim_avg_n_Wrap(data({year1:year2},:,:),0)
  rc2 = dim_avg_n_Wrap(data({year3:year4},:,:),0) - dim_avg_n_Wrap(data({year1:year2},:,:),0)
  return rc2
end 

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

function col_avg(data)
begin 

  plev = data&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)
  MSE_col =  dim_sum_wgt_n_Wrap(data, dz, 1, 1)
  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)
  MSE_nep = MSE_col
  MSE_nep = (/MSE_col/10^6/)

return MSE_nep
end 


;;;******************************************************************************
;;
begin

  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360

  yearstart = 198001
  yearend =202012

;;;;read data
    var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)     
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    q= q/1000.0

    hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

;;;

;;;;

    Cp = 1004
    Lv = 2.5*10^6
    g = 9.8
    
    CPT = t
    CPT = Cp*t
    
    lvq = q
    lvq = Lv*q
    
    gz = hgt
    gz = g*hgt

    MSE = t
    MSE = (/Cp*t + Lv*q + g*hgt/)
 
;;;MSE 
    MSE_nep = col_avg(MSE)
    
    CPT_nep = col_avg(CPT)

    lvq_nep = col_avg(lvq)

    gz_nep = col_avg(gz)
   


    MSE_nep = (/rmMonAnnCycTLL(MSE_nep)/) 
    CPT_nep = (/rmMonAnnCycTLL(CPT_nep)/)
    lvq_nep = (/rmMonAnnCycTLL(lvq_nep)/)
    gz_nep = (/rmMonAnnCycTLL(gz_nep)/)

;;;;;

;;;;;;split to 12 month  is a array size (12,yearnum)
   MSE_nep_12 = split_12(MSE_nep)
   CPT_nep_12 = split_12(CPT_nep)
   lvq_nep_12 = split_12(lvq_nep)
   gz_nep_12 = split_12(gz_nep)   
;;;;;


;;;等到年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)

    MSE_nep_annual_nep = dim_avg_n_Wrap(MSE_nep_annual(:,{0:20},{250:280}), (/1,2/))

   rc1 =regCoef_n(ispan(1,41,1),MSE_nep_annual_nep,0,0)

   rc1 = rc1*41
   diff  = (/avg(MSE_nep_annual_nep(20:))-avg(MSE_nep_annual_nep(0:19))/)
   print(diff)
exit()

;;;;;


;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta = (/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
    delete(date)
    delete(timeind)

    sst_ct = ssta
    sst_ct = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)


    sst_ct_djf = month_to_season(sst_ct, "DJF")
    ; sst_ct_djf(0:dimsizes(sst_ct_djf&time)-2,:,:) = sst_ct_djf(1:dimsizes(sst_ct_djf&time)-1,:,:)
    ; sst_ct_djf(dimsizes(sst_ct_djf&time)-1,:,:) = sst_ct(dimsizes(sst_ct&time)-1,:,:)

    nino3_4_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{190:240}), (/1,2/))   
    nino4_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{150:200}), (/1,2/))   
    nino3_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{210:270}), (/1,2/))   

    nino4_ct_sd = dim_standardize_n(nino4_ct, 1, 0)
    nino3_ct_sd = dim_standardize_n(nino3_ct, 1, 0)

    std_ct_nino3 = dim_stddev_n(nino3_ct, 0)
    std_ct_nino3_4 = dim_stddev_n(nino3_4_ct, 0)
   
    std_ct_nino4 = dim_stddev_n(nino4_ct, 0)
    
    ;;;;;
    

    sst_year = ispan(1980,2020,1)

    enino_year_ct_or = sst_year(ind(nino3_4_ct.ge.1.0*std_ct_nino3_4))
    lanina_year_ct_or = sst_year(ind(nino3_4_ct.le.-1.0*std_ct_nino3_4))
    
    enino_lanina_year = sst_year(ind(nino3_4_ct.ge.1.0*std_ct_nino3_4.or.nino3_4_ct.le.-1.0*std_ct_nino3_4))

    print(enino_year_ct_or)
    print(lanina_year_ct_or)
  
    ; ;;;for lanina
    ; enino_year_1 =(/1989,1999/)-1
    ; enino_year_2 =(/2000,2008,2011,2012/) -1

    ; enino_year_1 = enino_year_ct_or(0:3)
    ; lanina_year_1 = lanina_year_ct_or(0:2)

    enino_year_1 = enino_lanina_year(0:6)
    lanina_year_1 = enino_lanina_year(7:)


    enino_year1 = get1Dindex(sst_year,enino_year_1)
    lanina_year1 = get1Dindex(sst_year,lanina_year_1)




    ; ely=[1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016];
    ; enino_year = [/enino_year_ct,enino_year_wnp/]
;;;

;;;composite mean 
   MSE_enino = dim_avg_n_Wrap(MSE_nep_annual(enino_year1,:,:),0)
   MSE_lanina = dim_avg_n_Wrap(MSE_nep_annual(lanina_year1,:,:),0)
;;;;
   
   nino3_4_ct_3d = conform_dims(dimsizes(MSE_nep_annual), nino3_4_ct, (/0/))
   MSE_nep_annual_nor = MSE_nep_annual
   MSE_nep_annual_nor = (/MSE_nep_annual/abs(nino3_4_ct_3d)/)

   MSE_enino_nor = dim_avg_n_Wrap(MSE_nep_annual_nor(enino_year1,:,:),0)
   MSE_lanina_nor = dim_avg_n_Wrap(MSE_nep_annual_nor(lanina_year1,:,:),0)


   ; MSE_asy = MSE_enino
   ; MSE_asy = (/MSE_enino + MSE_lanina/)

   ; MSE_asy_nor = MSE_enino
   ; MSE_asy_nor = (/MSE_enino_nor - MSE_lanina_nor/)
   
;;;;

;;;   
  printMinMax(MSE_enino, 1)
  printMinMax(MSE_lanina, 1)
;;;
  
  MSE_DIFF = MSE_enino
  MSE_DIFF = (/MSE_lanina - MSE_enino/)
printVarSummary(MSE_DIFF)
  change = avg(MSE_DIFF({0:20},{250:280}))

  print(rc1)
  print(change)
  PERCENT= change/rc1
  print(PERCENT)
  exit()


  
  ; wks   = gsn_open_wks("pdf", "all-normalized-ENSO-ASY-MSE-spatial-pattern")
  ; wks   = gsn_open_wks("pdf", "all-elnino+lanina-ENSO-ASY-MSE-spatial-pattern")
  wks   = gsn_open_wks("pdf", "two-epoch-all-elnino+lanina-ENSO-ASY-MSE-spatial-pattern")

;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=0.015
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.2

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn= True
  ; res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360


  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.
  res@lbAutoManage = False
  res@lbLabelFontHeightF = 0.3

  res@tmXBValues = (/0,60,120,180,240,300,360/) 
  res@tmXBLabels = (/"0","60E","120E","180","120W","60W","0"/) 

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015

  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2

  res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  res@cnLevels = ispan(-10, 10, 2)*0.05


  res@lbLabelFontHeightF = 0.01
  res@lbLabelStride = 1
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueDarkRed18"
  res@cnFillPalette = "BlueWhiteOrangeRed"

  res@mpFillOn = True
  res@mpFillColor = "White"
  res@mpFillDrawOrder ="PostDraw"

  res@gsnRightString = ""

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpWidthF = 0.6
  ; res@vpHeightF = 0.1

    ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(4,graphic)

  ; res@tiMainString = data+ " MSE trend " + startyear + "-" + endyear
  res@gsnLeftString="JFMAM MSE&PRE El Nino + La Nina (1980-1999)" 
  plot(0) = gsn_csm_contour_map(wks, MSE_enino, res)
  res@gsnLeftString="JFMAM MSE&PRE El Nino + La Nina (2000-2020)" 
  plot(1) = gsn_csm_contour_map(wks, MSE_lanina, res)
  ; res@cnLevels = ispan(-10, 10, 2)*0.05
  res@cnLevels = ispan(-10, 10, 2)*0.05
; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)*0.25
  res@gsnLeftString="JFMAM MSE&PRE El Nino + La Nina (normalized) (1980-1999)" 
  plot(2) = gsn_csm_contour_map(wks, MSE_enino_nor, res)
  res@gsnLeftString="JFMAM MSE&PRE El Nino + La Nina (normalized) (2000-2020)" 
  plot(3) = gsn_csm_contour_map(wks, MSE_lanina_nor, res)

 
  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0

  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;
  add = "str" + ispan(1,4,1)
  
  do i=0,3
  plot@$add(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 

  ;;;打点
  ; ; end do
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnFillOn = True

  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLevels = (/0,0.1/)
  res_pro@lbLabelBarOn = False
  res_pro@cnMonoFillColor = False
  res_pro@cnFillColors = (/"Transparent","green","Transparent"/)


  res_pro@cnLinesOn=False
  ; res_pro@cnMonoLineColor =False
  ; res_pro@cnLineColors = ((/"Transparent","black","Transparent"/))

  res_pro@cnMonoFillPattern = True
  res_pro@cnFillPattern = 17
  res_pro@cnFillDotSizeF = 0.002
  
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
  ; iplot = new(2,graphic)
  ; iplot(0) = gsn_csm_contour(wks,prob_enino,res_pro)
  ; iplot(1) = gsn_csm_contour(wks,prob_lanina,res_pro)
  
;   do i=0,1 
;   overlay(plot(i), iplot(i))
; end do

    ;;  
   
pres = True
pres@gsnPanelLabelBar  = False
pres@pmLabelBarOrthogonalPosF = -0.01
gsn_panel(wks, plot, (/2,2/), pres)

end


