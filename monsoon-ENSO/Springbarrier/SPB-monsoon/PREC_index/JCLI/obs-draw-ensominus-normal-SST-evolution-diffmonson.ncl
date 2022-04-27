;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var = dim_variance_n(pv_ano(east_time,:,:), 0)
    var  =where(var.eq.0, var@_FillValue, var)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
    var@dim = n
 return var
 end

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
function composite_mean(ssta_month_12,u_850_month_12,v_850_month_12,ind1,ind2,ind3)
begin 

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))

    prob_ssta_str = one_test(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)
    prob_u_str = one_test(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    prob_v_str = one_test(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,u_month_str,v_month_str,prob_ssta_str,prob_u_str,prob_v_str/]
end 


;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 187101
    endmon = 201912
    startyear = 1871
    endyear = 2019

;;;读入HadiSST海温数据
    diri="/home/ys17-23/Extension/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
     
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_jja = nino34_djf
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

    nino34_djf = nino34_jja

    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear, 1)

;;;


;;;读入ERA wind data

;     ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
;     fils1 = "/home/ys17-23/Extension/20thCRV2/uwnd.mon.mean.nc"
;     f1 = addfile(fils1, "r")
;     time1 = f1->time
;     date1=cd_calendar(time1, 1)

;     endmon = endmon                 ;;; 1871-2012
;     startmon = startmon 
;     timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

 
; ;;;;U,V
;     u_200  = f1->uwnd(timeind1,{200},:,:)
;     u_850  = f1->uwnd(timeind1,{850},:,:)

;     u_850 = (/rmMonAnnCycTLL(u_850)/)
;     u_200 = (/rmMonAnnCycTLL(u_200)/)
; ;;;



;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("../All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2012})  ;;;1871-2012

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2012-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;


    ; f2=addfile("/home/ys17-23/Extension/GPCP/precip.mon.mean.nc", "r")
    ; date=cd_calendar(f2->time, 1)
    ; timeind=ind(date.ge.197901.and.date.le.201612)
    ; p=f2->precip(timeind,:,:)
   
    ; p = (/rmMonAnnCycTLL(p)/)
    ; p_ind = dim_avg_n_Wrap(p(:,{5:25},{60:100}), (/1,2/))
    ; p_ind_12 = split_12(p_ind)
    ; ; p_jjas = dim_avg_n_Wrap(p_ind_12(:,5:8), 1)
    ; std = dim_avg_n(dim_stddev_n(p_ind_12(:,5:8), 0),0)
    

    ; print(std)
    ; exit()

;;;strong and weak monsoon year 
    total_year = ispan(startyear+1,2011,1)   ;;;1872-2011
    DMI_ctrl_JJA_2 = DMI_ctrl_JJA(1:dimsizes(DMI_ctrl_JJA)-2)  ;;;1871-2011
    std = dim_stddev_Wrap(DMI_ctrl_JJA_2)


    inds = ind(DMI_ctrl_JJA_2.gt.0.43*std)
    indw = ind(DMI_ctrl_JJA_2.lt.-0.43*std)
    indnor = ind(DMI_ctrl_JJA_2.gt.-0.43*std.and.DMI_ctrl_JJA_2.lt.0.43*std)

    strmonrainfall = dim_avg_n_Wrap(DMI_ctrl_JJA_2(inds),0)/120.0
    print(strmonrainfall)
 

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal



 nino34_djf_2 = nino34_djf(1:dimsizes(nino34_djf)-1-(endyear-2011))


 wet_warm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;




;;;读入ERA wind data
    ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/ys17-23/Extension/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1871-2012
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  
;;;;
 
;;;;U,V

    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/ys17-23/Extension/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)
;;;;


;;;
    fils2 = "/home/ys17-23/Extension/20thCRV2/vwnd.mon.mean.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2 = cd_calendar(time2, 1)

    timeind2 = ind(date2.ge.startmon.and.date2.le.endmon)  ;;;风场数据从1979年开始

;;;;U,V
    v_200  = f2->vwnd(timeind1,{200},:,:)
    v_850  = f2->vwnd(timeind1,{850},:,:)

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)
   
    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
;;;

    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)
    dim_u = dimsizes(u_850)


;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


    case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

    case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

    case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]



    case = [/case1,case2,case3/]


  wksstring = (/"ElNino","LaNina","Normal"/)


do k=0,2 

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])
  
   u_year = ispan(startyear,endyear,1) 

   ind_us1_or :=  get1Dindex(u_year,case_need[0])
   ind_uw1_or :=  get1Dindex(u_year,case_need[1])
   ind_unor1 :=  get1Dindex(u_year,case_need[2])
   

   ind_uw1 :=  ind_uw1_or(0:dimsizes(ind_uw1_or)-2)  ;;去掉最后一年

   ind_us1 :=  ind_us1_or(0:dimsizes(ind_us1_or)-2)  ;;去掉最后一年


  plot = new(6,graphic)
  plot1 = new(6,graphic)
  plot2 = new(6,graphic)
  iplot = plot1
  iplot1 = plot
  iplot2 = plot


  ; wks = gsn_open_wks("pdf", "seasonalmean-WNP-UVSST_strong-monsoon-after")
  ; wks1 = gsn_open_wks("pdf", "seasonalmean-WNP-UVSST_weak-monsoon-after")

  wks = gsn_open_wks("pdf", "diff-9month-seasonalmean-UVSST-Str-normon-"+" "+wksstring(k))
  wks1 = gsn_open_wks("pdf", "diff-9month-seasonalmean-UVSST-Weak-normon-"+" "+wksstring(k))
  

 ; month = (/"JJA","SON","DJF(+1)","MAM(+1)"/)  
  month = (/"MAM","JJA","SON","DJF(+1)","MAM(+1)","JJA(+1)"/) 
 
 
;;;split to 12 month 
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 -1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  dim2  = dimsizes(u_850)
  yearnum := dim2(0)/12 -1
  u_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  v_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)


 do i=0,11   
     if i.le.9 then
       ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1-12:12,:,:)
       u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1-12:12,:,:)
       v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1-12:12,:,:)
     else
        ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1:12,:,:)
        u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1:12,:,:)
        v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1:12,:,:)
     end if
 end do 


do i=0,5
;;; strong and weak monsoon year ssta , u ,v 
   if i.le.3 then 
   
    com_month_str  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_ts1,ind_us1,ind_us1)
    
    ssta_month_str = com_month_str[0]
    u_month_str = com_month_str[1]
    v_month_str = com_month_str[2]

    prob_ssta_str = com_month_str[3]
    prob_u_str = com_month_str[4]
    prob_v_str = com_month_str[5]

    com_month_weak  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_tw1,ind_uw1,ind_uw1)
    
    ssta_month_weak = com_month_weak[0]
    u_month_weak = com_month_weak[1]
    v_month_weak = com_month_weak[2]

    prob_ssta_weak = com_month_weak[3]
    prob_u_weak = com_month_weak[4]
    prob_v_weak = com_month_weak[5]

    com_month_nor  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_tnor1,ind_unor1,ind_unor1)
    
    ssta_month_nor = com_month_nor[0]
    u_month_nor = com_month_nor[1]
    v_month_nor = com_month_nor[2]

    prob_ssta_nor = com_month_nor[3]
    prob_u_nor = com_month_nor[4]
    prob_v_nor = com_month_nor[5]
    else 
      
    com_month_str  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_ts1+1,ind_us1+1,ind_us1+1)
    
    ssta_month_str = com_month_str[0]
    u_month_str = com_month_str[1]
    v_month_str = com_month_str[2]

    prob_ssta_str = com_month_str[3]
    prob_u_str = com_month_str[4]
    prob_v_str = com_month_str[5]

    com_month_weak  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_tw1+1,ind_uw1+1,ind_uw1+1)
    
    ssta_month_weak = com_month_weak[0]
    u_month_weak = com_month_weak[1]
    v_month_weak = com_month_weak[2]

    prob_ssta_weak = com_month_weak[3]
    prob_u_weak = com_month_weak[4]
    prob_v_weak = com_month_weak[5]

    com_month_nor  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_tnor1+1,ind_unor1+1,ind_unor1+1)
    
    ssta_month_nor = com_month_nor[0]
    u_month_nor = com_month_nor[1]
    v_month_nor = com_month_nor[2]

    prob_ssta_nor = com_month_nor[3]
    prob_u_nor = com_month_nor[4]
    prob_v_nor = com_month_nor[5]
    
    end if

    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
 
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
; 
  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=False
  res@cnLinesOn=True

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 300

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@gsnContourNegLineDashPattern = 1
  res@gsnContourZeroLineThicknessF = 2.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.2

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  
  res@gsnLeftString = month(i) + " SST&850UV " 
  res@gsnLeftStringFontHeightF = 0.015

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

   ; res@gsnScalarContour        = True

    ; ;;;;;;;LineArrow
    ; res@vcMinDistanceF           = 0.02           ; thin vectors
    ; res@vcMinMagnitudeF          =0.4
    ; res@vcRefMagnitudeF          = 4.0
    ; res@vcRefLengthF             = 0.1     ; define length of v
    ; ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ; res@vcLineArrowThicknessF = 1.0
    ; ; res@vcLineArrowHeadMaxSizeF=0.008
    ; ; res@vcLineArrowHeadMinSizeF=0.008
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006

    ; if i.eq.5 then 
    ; res@vcRefAnnoOn = True
    ; res@vcRefAnnoOrthogonalPosF = 0.2
    ; else 
    ; res@vcRefAnnoOn = False
    ; end if

    u_month_str = (/u_month_str - u_month_nor/)
    v_month_str = (/v_month_str - v_month_nor/)
    ssta_month_str = (/ssta_month_str - ssta_month_nor/)

    u_month_weak = (/u_month_weak - u_month_nor/)
    v_month_weak = (/v_month_weak - v_month_nor/)
    ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)

    var1 = (/prob_ssta_str/)
    var2 = (/prob_ssta_nor/)
    var3 = (/prob_ssta_weak/)



    s1 = prob_ssta_str@dim
    s2 = prob_ssta_nor@dim
    s3 = prob_ssta_weak@dim
  

    prob_ssta_str = (/ttest(ssta_month_str, var1,s1 , ssta_month_nor, var2, s2, False, False)/)
    prob_ssta_weak = (/ttest(ssta_month_weak, var3, s3, ssta_month_nor, var2, s2, False, False)/)
    
    copy_VarCoords(ssta_month_str,prob_ssta_str)
    copy_VarCoords(ssta_month_str,prob_ssta_weak)


    res@gsnRightString = wksstring(k) + " event (Str Mon - Without Mon)"
    ; res@gsnRightString = "Strong ISM"
    ; plot(i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
    plot(i) = gsn_csm_contour_map(wks,ssta_month_str, res)
    res@gsnRightString = wksstring(k) + " event (Weak Mon - Without Mon)"
    ; res@gsnRightString = "Weak ISM"

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -0.5
    ; res@cnMaxLevelValF       =   0.5
    ; res@cnLevelSpacingF      =  0.1

    ; plot1(i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)
    plot1(i) = gsn_csm_contour_map(wks,ssta_month_weak, res)
  
  
  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0
  ;;;画box
  lat1 = -5
  lat2 = 5
  lon1 = 190
  lon2 = 240

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;

   plot@$month(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   plot1@$month(i)$ = gsn_add_polyline(wks1, plot1(i), pgx_pool, pgy_pool, pgres)
   ; plot2@$month(i)$ = gsn_add_polyline(wks2, plot2(i), pgx_pool, pgy_pool, pgres)

   ; str1 = unique_string("dum1")  ; "unique_string" will return a unique    

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.05/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.001
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

  
    iplot(i) = gsn_csm_contour(wks,prob_ssta_str,res_pro)
    iplot1(i) = gsn_csm_contour(wks1,prob_ssta_weak,res_pro)
    ; iplot2(i) = gsn_csm_contour(wks2,prob_ssta_nor,res_pro)

    ; overlay(plot(i), iplot(i))
    ; overlay(plot1(i), iplot1(i))
    ; overlay(plot2(i), iplot2(i))
    delete(prob_ssta_str)
    delete(prob_ssta_weak)
end do
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelStride=1
; gsn_panel(wks, plot, (/3,2/), pres)
; gsn_panel(wks1, plot1, (/3,2/), pres)
; gsn_panel(wks2, plot2, (/3,2/), pres)
gsn_panel(wks, plot(0:3), (/2,2/), pres)
gsn_panel(wks1, plot1(0:3), (/2,2/), pres)
end do
end 