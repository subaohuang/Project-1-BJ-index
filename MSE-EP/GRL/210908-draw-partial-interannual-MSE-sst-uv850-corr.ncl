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

function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

function part_corr(var1,var2,var3)
begin 
  ;;;remove var3 influence 
  r12 = escorc_n(var1, var2, 0, 0)
  ; copy_VarCoords(var2(0,:,:), r12)  ;;;这是二维
  r13 = escorc_n(var1, var3, 0, 0)  ;;;这个是1dim
  r23 = escorc_n(var2, var3, 0, 0)   ;;;这个是2dim
  r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
  rc_mam = r12
  rc_mam = (/r12_3/)
return rc_mam
end 

function band_pass(data)
begin
  ihp   = 2                            ; band_pass
  sigma = 1.0                           ; Lanczos sigma
  nWgt  = 9                           ; will lose nwgt-1 data must odd number 
  fcb   = 1./2.                        ; 2 years
  fca   = 1./9.                         ;7年
  wgts  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
  msf_500_band=wgt_runave_Wrap(data, wgts, 0)
  return msf_500_band
  end 
begin
    ;;;time span
    startmon =198001
    endmon = 202012

    running_time = 15

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;; 
;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon-100)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta = sst
    ssta=(/rmMonAnnCycTLL(sst)/)


    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    copy_VarCoords(sst, ssta)
   ;;
;;;
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_ncep1_50-20.nc", "r")

    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    startmon = 198001
    endmon = 202012

    startyear = 1980
    endyear = 2020

    timeind :=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)
   

    ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)
   

    ; MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_nep = MSE_col
    MSE_nep = (/MSE_col/10^6/)
;;;;


;;;;
  path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
  var_name = (/"uwnd.mon.mean.nc","vwnd.mon.mean.nc","shum.mon.mean.nc","air.mon.mean.nc"/)

  filein = path_wind + var_name 

  ;;;region 
  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360
  ;;;

  ; ;;;global region 
  ; lat1 = -50
  ; lat2 = 50
  ; lon1 = 0
  ; lon2 = 360
  ; ;;;
  startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
  endyear_r = str_split_by_length(endmon, 4 )

  startyear = stringtoint(startyear_r(0))
  endyear = stringtoint(endyear_r(0))

  yearnum = endyear - startyear +1


  f1 = addfile(filein(0),"r")

  date:=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)


  u = f1->uwnd(timeind,{850},{lat1:lat2},{lon1:lon2})
  
  u = (/rmMonAnnCycTLL(u)/)
  u = (/dtrend_msg_n(ispan(1,dimsizes(u&time),1),u,False,True,0)/)

  f2 = addfile(filein(1),"r")
  date:=cd_calendar(f2->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  v = f2->vwnd(timeind,{850},{lat1:lat2},{lon1:lon2})

  v = (/rmMonAnnCycTLL(v)/)
  v = (/dtrend_msg_n(ispan(1,dimsizes(v&time),1),v,False,True,0)/)
;;;;
  
  f3 = addfile(filein(2),"r")
  ; f3 = addfile(filein(3),"r")
  date:=cd_calendar(f3->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)
  q = f3->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

  ; f4 = addfile(filein(3),"r")
  ; f3 = addfile(filein(3),"r")
  ; date:=cd_calendar(f4->time, 1)
  ; timeind:=ind(date.ge.startmon.and.date.le.endmon)
  ; t = f4->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

  Cp = 1004
  Lv = 2.5*10^6
  g = 9.8

  ; t = t +273.15
  ; t = Cp*t 
  ; q = (/t + Lv*q/)

  plev := q&level
  ptop = min(plev)
  psfc = 1013.0
  dz := dpres_plevel(plev, psfc, ptop, 0)

  q_col =  dim_sum_wgt_n_Wrap(q, dz, 1, 1)

  q_col = (/rmMonAnnCycTLL(q_col)/)

  q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
   
  ; MSE_nep = q_col


;;;;;;split to 12 month  is a array size (12,yearnum)

   MSE_nep_12 = split_12(MSE_nep)
   ssta_12 = split_12(ssta)

   u_12 = split_12(u)
   v_12 = split_12(v)
;;;;;

   
;;;等到年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual!0 = "time"
    MSE_nep_annual&time = ispan(startyear,endyear,1)

    ssta_annual = month_to_annual(ssta, 1)
    ssta_annual!0 = "time"
    ssta_annual&time = ispan(startyear-1,endyear,1)
    ssta_annual = (/dim_avg_n(ssta_12(0:4,:,:,:),0)/)

    ssta_mam = ssta_annual(1:dimsizes(ssta_annual&time)-1,:,:)


    
   ; ;;;;去掉前期DJF的SST信号
   ;  ssta_annual = month_to_season(ssta, "DJF")
   ;  ssta_annual&time = ispan(startyear-1,endyear,1)
   ;  ssta_mam = ssta_annual(1:dimsizes(ssta_annual&time)-1,:,:)
   ; ;;;;



    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)

    u_annual = month_to_annual(u, 1)
    v_annual = month_to_annual(v, 1)

    u_annual = (/dim_avg_n(u_12(0:4,:,:,:),0)/)
    v_annual = (/dim_avg_n(v_12(0:4,:,:,:),0)/)
;;;;;

;;;
    MSE_nep_mam = dim_avg_n_Wrap(MSE_nep_annual(:,{0:20},{250:280}), (/1,2/))

    ;;;局地海温
    ; ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

    ;  ;;nino1+2:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{270:280}), (/1,2/))
    ; ;;;
    
    ;  ;;nino3:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{150:270}), (/1,2/))
    ; ;;;
     
     ;;;nino4 
      nino3_4 = dim_avg_n_Wrap(ssta_mam(:,{-5:5},{190:240}), (/1,2/))
      ; nino3_4_decadal = runave_n_Wrap(nino3_4, 11, 0, 0)
      ; nino3_4 = (/nino3_4 - nino3_4_decadal/)
      ; nino3_4 = (/band_pass(nino3_4)/)
      atlantic = dim_avg_n_Wrap(ssta_mam(:,{10:30},{280:330}), (/1,2/))
      ; atlantic = dim_avg_n_Wrap(ssta_mam(:,{0:20},{40:100}), (/1,2/))  ;;;indian
      ; atlantic = dim_avg_n_Wrap(ssta_mam(:,{-20:0},{300:360}), (/1,2/))  ;;;indian
      ; atlantic = (/band_pass(atlantic)/)
    ;;;

    

    ;;;;去掉局地海温的影响,偏相关性
    ;;; 1 for MSE, 2 for global sst , 3 for nino3_4 ,4 for atlantic SST
        r12no3 = part_corr(MSE_nep_mam,ssta_mam, nino3_4) ;;;12.3  2dim
        ; rc_mam_rmaltantic = part_corr(MSE_nep_mam,ssta_mam,nino3_4) ;;;12.4 
        r12no4 = part_corr(MSE_nep_mam,ssta_mam, atlantic) ;;; 2dim

        r13no4 = part_corr(MSE_nep_mam,nino3_4, atlantic)  ;;; 1dim

        r23no4 = part_corr(ssta_mam,nino3_4,atlantic)    ;;; 1dim
       
        rc_mam = (/r12no3/)

        r13no4_23no4 = r13no4*r23no4

        r13no4_23no4_2dim  = r12no4
        r13no4_23no4_2dim  = r13no4_23no4
        ; rc_mam = (/r12no4 - r13no4_23no4_2dim/)/(sqrt((1- r13no4^2)*(1- r23no4^2)))

        copy_VarCoords(ssta_mam(0,:,:), rc_mam)
    ;;;;
    


    ;;;;去掉局地海温的影响,偏相关性
    ;;; 1 for MSE, 2 for global u , 3 for nino3_4 ,4 for atlantic SST
        r12no3_u = part_corr(MSE_nep_mam, u_annual, nino3_4) ;;;12.3  2dim
        ; rc_mam_rmaltantic = part_corr(MSE_nep_mam,ssta_mam,nino3_4) ;;;12.4 
        r12no4_u = part_corr(MSE_nep_mam, u_annual, atlantic) ;;; 2dim

        r13no4_u = part_corr(MSE_nep_mam, nino3_4, atlantic)  ;;; 1dim

        r23no4_u = part_corr(u_annual, nino3_4,atlantic)    ;;; 1dim
       
        rc_mam_u = (/r12no3_u/)

        r13no4_23no4_u = r13no4_u*r23no4_u

        r13no4_23no4_u_2dim  = r12no4_u
        r13no4_23no4_u_2dim  = r13no4_23no4_u
        ; rc_mam_u = (/r12no4_u - r13no4_23no4_u_2dim/)/(sqrt((1- r13no4_u^2)*(1- r23no4_u^2)))

        copy_VarCoords(u_annual(0,:,:), rc_mam_u)
    ;;;;
    

    ;;;;
        ;;;;去掉局地海温的影响,偏相关性
    ;;; 1 for MSE, 2 for global u , 3 for nino3_4 ,4 for atlantic SST
        r12no3_v = part_corr(MSE_nep_mam, v_annual, nino3_4) ;;;12.3  2dim
        ; rc_mam_rmaltantic = part_corr(MSE_nep_mam,ssta_mam,nino3_4) ;;;12.4 
        r12no4_v = part_corr(MSE_nep_mam, v_annual, atlantic) ;;; 2dim

        r13no4_v = part_corr(MSE_nep_mam, nino3_4, atlantic)  ;;; 1dim

        r23no4_v = part_corr(v_annual, nino3_4,atlantic)    ;;; 1dim
       
        rc_mam_v = (/r12no3_v/)

        r13no4_23no4_v = r13no4_v*r23no4_v

        r13no4_23no4_v_2dim  = r12no4_v
        r13no4_23no4_v_2dim  = r13no4_23no4_v
        ; rc_mam_v = (/r12no4_v - r13no4_23no4_v_2dim/)/(sqrt((1- r13no4_v^2)*(1- r23no4_v^2)))

        copy_VarCoords(v_annual(0,:,:), rc_mam_v)
 ;;;;;;


    ;;;;有效自由度
      
        N = dimsizes(ssta_mam(:,1,1))

        ; effN = N-2  ;;; for remove two variable 
        effN = N-1
        ; prob_mam = rtest(rc_mam, N-1, 0)  ;;;1个变量去掉，加1
        prob_mam = rtest(rc_mam, effN, 0)  ;;;1个变量去掉，加1
        copy_VarCoords(ssta_mam(0,:,:), prob_mam)

        prob_mam_u = rtest(rc_mam_u, effN, 0)  ;;;1个变量去掉，加1
        copy_VarCoords(u_annual(0,:,:), prob_mam_u)

        prob_mam_v = rtest(rc_mam_v, effN, 0)  ;;;1个变量去掉，加1
        copy_VarCoords(v_annual(0,:,:), prob_mam_v)

    ;;;;;
    


    ; ;;;normal corr

    ;   rc_mam = escorc_n(MSE_nep_mam, ssta_mam, 0, 0)

    ;   copy_VarCoords(ssta_mam(1,:,:), rc_mam)

    ;   prob_mam = rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
    ;   copy_VarCoords(ssta_mam(0,:,:), prob_mam)

    ;   rc_mam_u = escorc_n(MSE_nep_mam, u_annual, 0, 0)

    ;   copy_VarCoords(u_annual(1,:,:), rc_mam_u)

    ;   prob_mam_u = rtest(rc_mam_u, dimsizes(ssta_mam&time)-1, 0)
    ;   copy_VarCoords(u_annual(0,:,:), prob_mam)


    ;   rc_mam_v = escorc_n(MSE_nep_mam, v_annual, 0, 0)

    ;   copy_VarCoords(v_annual(1,:,:), rc_mam_v)

    ;   prob_mam_v = rtest(rc_mam_v, dimsizes(ssta_mam&time)-1, 0)
    ;   copy_VarCoords(v_annual(0,:,:), prob_mam_v)

    ;   ;;;;


    plot = new(4,graphic)
;;;;

    rc_mam_u = where(prob_mam_u.le.0.05.or.prob_mam_v.le.0.05, rc_mam_u, rc_mam_u@_FillValue)
    rc_mam_v = where(prob_mam_u.le.0.05.or.prob_mam_v.le.0.05, rc_mam_v, rc_mam_u@_FillValue)
    
; wks   = gsn_open_wks("pdf", "MSE与SST的相关-JFMAM-remove_Nino3.4_rmi(ndian")
; wks   = gsn_open_wks("pdf", "MSE与SST&UV850的相关-JFMAM-remove_Nino3.4 (Interannual)")
; wks   = gsn_open_wks("pdf", "Lvq与SST&UV850的相关-JFMAM-remove_Nino3.4 & NTA (Interannual)")
   ; wks   = gsn_open_wks("pdf", "去掉前期ENSOLvq与SST&UV850的相关-JFMAM-remove_Nino3.4 & NTA")
   wks   = gsn_open_wks("pdf", "去掉前期ENSO-MSE与SST&UV850的相关-JFMAM-remove_Nino3.4")
; wks   = gsn_open_wks("pdf", "MSE与SST&UV850的相关-JFMAM-remove_Nino3.4_rmatlantic")
; wks   = gsn_open_wks("pdf", "MSE与SST&UV850的相关-JFMAM-remove_Nino3.4_rmatlantic")

; wks   = gsn_open_wks("pdf", "Interannual-MSE与SST&UV的相关-JFMAM")
; wks   = gsn_open_wks("pdf", "q与SST&UV的相关-JFMAM")
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


  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2
  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  res@lbLabelStride = 1

  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"

  res@gsnScalarContour        = True
 


  ; res@vcGlyphStyle="CurlyVector"
  res@vcGlyphStyle="LineArrow"
  res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.015          ; thin vectors

  res@vcRefLengthF             = 0.05    ; define length of v
  res@vcVectorDrawOrder="PostDraw"

  ; resv@vcMonoLineArrowColor = True
  ; resv@vcLineArrowColor = "blue"
  res@vcLineArrowThicknessF=1.2
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  res@vcRefMagnitudeF          = 0.2*10
  ; res@vcMinMagnitudeF          =.3

  res@vcRefAnnoString1On = True
  res@vcRefAnnoString1 = "2.0 m/s"
  res@vcRefAnnoOrthogonalPosF  = -1.21              ; move ref vector down -1.0在图的右上
  res@vcRefAnnoParallelPosF = 1.02
  res@vcRefAnnoFontHeightF = 0.013
  res@vcRefAnnoBackgroundColor = "Transparent"
  res@vcRefAnnoPerimOn = False


   ; res@gsnLeftString="Corr (NEP MSE & SST UV850)"  
   ; res@gsnLeftString="Corr (NEP MSE & SST UV850)"  
  ; res@gsnLeftString="Partial corr (NEP MSE & SST UV850, remove Nino3.4 & IND SST)" 
   ; res@gsnLeftString="Partial corr (NEP MSE & SST UV850, remove Nino3.4 & NTA SST)" 
  ; res@gsnLeftString="Partial corr (NEP  CpT & SST UV850, remove Nino3.4 & NTA (Interannual))" 
  res@gsnLeftString="Partial corr (NEP  CpT & SST UV850, remove Nino3.4 & NTA )" 
  ; res@gsnLeftString="Partial corr (NEP MSE & SST , remove Nino3.4 & NTA SST)"

  res@gsnLeftStringFontHeightF = 0.015
  ; res@gsnLeftString = "The components of JFMAM shum & 850UV trend linked to TNA"
  res@gsnLeftStringOrthogonalPosF  = -0.13
 
  plot = new(4,graphic)

  plot(0) = gsn_csm_vector_scalar_map(wks,rc_mam_u,rc_mam_v,rc_mam, res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    



  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=2.0
  pgres@gsLineDashPattern=0


  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;

  str2 = ispan(21, 24, 1) +""
  do i =0,3 
   plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(4,graphic)
     i=0
    iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)
  
  
   
     overlay(plot(i), iplot(i))
  ; end do

    ;;;  
   
; pres = True
; pres@gsnPanelLabelBar  = True
; pres@pmLabelBarOrthogonalPosF = -0.01
draw(plot(0))
frame(wks)
; gsn_panel(wks, plot, (/2,2/), pres)

end


