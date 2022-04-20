; function  rc_test(ssta_rc)
; begin
;   tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
;    df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
;    b = tval    ; b must be same size as tval (and df)
;    b = 0.5
;    prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
;    copy_VarCoords(ssta_rc, prob)

; return prob
; end
;;;add ncl_to_py
function  rc_test(ssta_rc,id1,id2)
begin
   dim = dimsizes(id1)
   varname = getVarDimNames(id2)
   ;;有效自由度 
   r1_r = esacr_n(id1, 1,0)  
   r2_r = esacr_n(id2, 1,0)

   r1_o = r1_r(1)
   r2 = r2_r(1,:,:)  ;;;取lag1
   f = dimsizes(r2)

   r1  = onedtond(r1_o, dimsizes(r2))

   v = dim*(/1-r1*r2/)/(/1+r1*r2/) 
 
   ssta_rc_new = onedtond(ssta_rc, dimsizes(ssta_rc))
   tval = sqrt(v)*(/ssta_rc_new/)/(sqrt(1- ssta_rc_new^2))

   df   = v
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
   ;;或者用student_t
  ;; prob_1 = student_t(tval, df)x
return prob
end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

begin
    yearstart = 194901
    yearend = 201912

    lat1 = -30
    lat2 = 30
    lon1 = 0
    lon2 = 360
    
    runmean_time = 11


    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
    MSE = f1->MSE

    ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

    ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

    ; MSE_col = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_col = (/MSE_col/10^6/)
    MSE_col_annual = month_to_annual(MSE_col, 1)
    MSE_col_annual!0 = "time"

    MSE_col_MAM = month_to_season(MSE_col, "MAM")
    MSE_col_JJA = month_to_season(MSE_col, "JJA")
    MSE_col_SON = month_to_season(MSE_col, "SON")
    MSE_col_DJF = month_to_season(MSE_col, "DJF")

    ;;;;;;split to 12 month nino34_12 is a array size (12,yearnum)
  
    dim = dimsizes(MSE_col) 
    yearnum = dim(0)/12

    MSE_col_12 = new((/12,yearnum,dim(1),dim(2)/),float)

    do i=0,11   
     MSE_col_12(i,:,:,:) = MSE_col(i::12,:,:)
    end do 

    MSE_col_MAM = (/dim_avg_n(MSE_col_12(0:4,:,:,:),0)/)

    MSE_nep_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{0:20},{240:280}), (/1,2/))
    
    MSE_nep_mam = (/runave_n(MSE_nep_mam, runmean_time, 0, 0)/)


    ;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201912.and.date.ge.195001)

    sst  = lonFlip(f1->sst(timeind,{-50:50},:))

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    ; ssta = (/runave_n_Wrap(ssta, 9*12, 0, 0)/)
    ssta_mam = month_to_season(ssta, "MAM")

       ;;;;;;split to 12 month nino34_12 is a array size (12,yearnum)
        
        dim := dimsizes(ssta)
        yearnum := dim(0)/12
        ssta_12 = new((/12,yearnum,dim(1),dim(2)/),float)

      do i=0,11   
       ssta_12(i,:,:,:) = ssta(i::12,:,:)
      end do 

    ssta_mam = (/dim_avg_n(ssta_12(0:4,:,:,:),0)/)

    ssta_mam = (/runave_n(ssta_mam, runmean_time, 0, 0)/)

    ;;;局地海温
    ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))
    ;;;
      
    ;;;;去掉局地海温的影响,偏相关性

    r13 = escorc_n(MSE_nep_mam, ssta_mam_nep, 0, 0)  ;;;这个是1dim
    r23 = escorc_n(ssta_mam_nep, ssta_mam, 0, 0)   ;;;这个是2dim

    ;;;计算相关系数
      rc_mam = escorc_n(MSE_nep_mam, ssta_mam, 0, 0)
      copy_VarCoords(ssta_mam(0,:,:), rc_mam)  ;;;这是二维
      printVarSummary(rc_mam)
      r12 = rc_mam
    ;;;

    ; r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
    
    ; rc_mam = (/r12_3/)

   ;;
    ; prob_mam = rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
    ; copy_VarCoords(ssta_mam(0,:,:), prob_mam)
    
     prob_mam = rc_test(rc_mam,MSE_nep_mam,ssta_mam)
    copy_VarCoords(ssta_mam(0,:,:), prob_mam)


    ;;偏相关的有效自由度
   ;   ; r1_auto = esacr_n(MSE_nep_mam,1,0)
   ;   r2_auto = esacr_n(ssta_mam,1,0)
   ;   exit()
   ;   one = r2_auto(0,:,:)
   ;   one =1
    
   ;   N = dimsizes(ssta_mam&time)*(one- r1_auto(0)*r2_auto(0,:,:))/(one+r1_auto(0)*r2_auto(0,:,:))
   ; print(N)
   ; exit()
    ;   prob_mam = rtest(rc_mam, N, 0)
    ; copy_VarCoords(ssta_mam(0,:,:), prob_mam)

    ; prob_mam = rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
    ; copy_VarCoords(ssta_mam(0,:,:), prob_mam)

    plot = new(4,graphic)
;;;;

    
; wks   = gsn_open_wks("pdf", "nodetrend-MSE与SST的相关-MAM-"+runmean_time+"年滑动平均")
wks   = gsn_open_wks("pdf", "MSE与SST的相关-JFMAM-"+runmean_time+"年滑动平均")
; wks   = gsn_open_wks("pdf", "MSE与SST的相关-MAM-去掉局地的影响")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=True

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


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"
  res@tiMainString = "Corr coef between" +runmean_time+"-year runmean SST and area-averaged NEP MSE during JFMAM"
  res@tiMainFontHeightF = 0.015

  res@gsnLeftString="" 
  ; res@gsnLeftString="(COR MAM MSE, MAM SST)" 
  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam ,res)
  ; res@gsnLeftString="(b)JJA MSE trend " 
  ; plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  ; res@gsnLeftString="(c)SON MSE trend " 
  ; plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  ; res@gsnLeftString="(d)DJF MSE trend " 
  ; plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 


   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 210
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
   do i =0,3 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 

    iplot=new(4,graphic)

  ;   ;;;打点1
  ;   res_pro=True
  ;   res_pro@gsnDraw=False
  ;   res_pro@gsnFrame=False
  ;   res_pro@cnLevelSelectionMode = "ExplicitLevels"
  ;   res_pro@cnLinesOn=False
  ;   res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ;   res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  ;   res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  ;   res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

  ;   res_pro@cnLineLabelsOn=False
  ;   res_pro@tmYRMode = "Automatic"
  ;   res_pro@cnInfoLabelOn=False


    
 
  ;   iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    
  ;   ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)    
  ;   ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)    
  ;   ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)    

  ;   ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
  ;   ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
  ;   ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

  ;   opt                      = True
  ;   opt@gsnShadeFillType     = "pattern"
  ;   opt@gsnShadeMid     = 17         ; stipple
  ;   ; opt@gsnShadeFillDotSizeF =  0.0025
  ;   ; opt@gsnShadeHigh   = "green"
  ;   opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

  ;   i=0
  ;    ; do i=0,3
  ;   iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
  ;    overlay(plot(i), iplot(i))
  ; ; end do
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnFillOn = True

  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLevels = (/0,0.1/)
  res_pro@lbLabelBarOn = False
  res_pro@cnMonoFillColor = False
  res_pro@cnFillColors = (/"Transparent","black","Transparent"/)


  ; res_pro@cnLinesOn=True
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
  i=0
  iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)

  overlay(plot(i), iplot(i))
  ; end do

; pres = True
; pres@gsnPanelLabelBar  = True
; pres@pmLabelBarOrthogonalPosF = -0.01
draw(plot(0))
frame(wks)
; gsn_panel(wks, plot, (/2,2/), pres)

end


