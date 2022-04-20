;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
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
function composite_mean(ssta_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},{lon1:lon2}),(/1,2,3/))

return [/ssta_month_str/]
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
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst=lonFlip(f1->sst(timeind,:,:))
   
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
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)  ;;;;JJA > 0.43*STD

    nino34_djf = nino34_jja

    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear, 1)

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))
;;;



;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear+1:2016})    ;;;1872-2012

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2016-startyear-1+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;

;;;strong and weak monsoon year 
    total_year = ispan(startyear+1,2016,1)

    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    inds = ind(DMI_ctrl_JJA.gt.0.43*std)
    indw = ind(DMI_ctrl_JJA.lt.-0.43*std)
    indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(1:dimsizes(nino34_djf)-1-(endyear-2016))   ;;;;1982-2016

 wet_warm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;


    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)
    ; dim_u = dimsizes(u_850)


;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"(a) El Nino","(b) La Nina","(c) Normal"/)


  plot = new(3,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "JFMAM2JFMAM-ssta-evolution-during-diffmonsoon")
   
  stss = ispan(1,3,1)+""
do k=0,2 
  print("***********start*********************"+k)

  case_need:= case[k]

  ind_ts1 :=  get1Dindex(sst_year,case_need[0])
  ind_tw1 :=  get1Dindex(sst_year,case_need[1])
  ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

  ;;;split to 12 month 
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)
  ;;;

  len  = 27

  ssta_str_evo = new((/len/),float)
  ssta_weak_evo = new((/len/),float)
  ssta_nor_evo = new((/len/),float)

 do i=0,11   
       ssta_month_12(i,:,:,:) = ssta(i:dimsizes(ssta&time)-1:12,:,:)   ;;;1月开始
 end do 

    com_month_str  = composite_mean(ssta_month_12(5:11,:,:,:),ind_ts1-1)
    
    ssta_str_evo(0:6) = com_month_str[0]

    com_month_weak  = composite_mean(ssta_month_12(5:11,:,:,:),ind_tw1-1)
    
    ssta_weak_evo(0:6) = com_month_weak[0]
   
    
    com_month_nor  = composite_mean(ssta_month_12(5:11,:,:,:),ind_tnor1-1)
   
    ssta_nor_evo(0:6) = com_month_nor[0]
 

    ;;;;
    com_month_str  = composite_mean(ssta_month_12,ind_ts1)
    
    ssta_str_evo(7:18) = com_month_str[0]
   


    com_month_weak  = composite_mean(ssta_month_12,ind_tw1)
    
    ssta_weak_evo(7:18) = com_month_weak[0]
  

    ; ; if sst_year(ind_tnor1(dimsizes(ind_tnor1)-1)).eq.141 then 
    ; ;   ind_tnor1_r = ind_tnor1(0:dimsizes(ind_tnor1)-2)
    ; ;   ind_unor1_r = ind_unor1(0:dimsizes(ind_unor1)-2)
    ; ; else 
    ;   ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)
    ;   ind_unor1_r := ind_unor1(0:dimsizes(ind_unor1)-2)
    ; ; end if
    
    com_month_nor  = composite_mean(ssta_month_12,ind_tnor1)
   
    ssta_nor_evo(7:18) = com_month_nor[0]
   
   
    com_month_str  = composite_mean(ssta_month_12(0:7,:,:,:),ind_ts1+1)
    
    ssta_str_evo(19:) = com_month_str[0]
   

    com_month_weak  = composite_mean(ssta_month_12(0:7,:,:,:),ind_tw1+1)
    
    ssta_weak_evo(19:) = com_month_weak[0]
    
    
    com_month_nor  = composite_mean(ssta_month_12(0:7,:,:,:),ind_tnor1+1)
    
    ssta_nor_evo(19:) = com_month_nor[0]
   

    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
 

;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  res@tmXTOn         = False
  res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = len - 1



  res@gsnYRefLine = 0.0


  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, len-1, 1)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
  res@tmXBLabels = (/"J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
   res@tmYLMode = "Explicit"



  if k.eq.0 then 
   res@trYMinF        = -0.6
   res@trYMaxF        = 1.6

   res@tmYLValues := ispan(-6, 16, 2)*0.1
   res@tmYLLabels := ispan(-6, 16, 2)*0.1
  end if 

    if k.eq.1 then 
   res@trYMinF        = -1.2
   res@trYMaxF        = 0.4

   res@tmYLValues := ispan(-12, 4, 2)*0.1
   res@tmYLLabels := ispan(-12, 4, 2)*0.1
  end if 

     if k.eq.2 then 
   res@trYMinF        = -0.6
   res@trYMaxF        = 0.6

   res@tmYLValues := ispan(-6, 6, 2)*0.1
   res@tmYLLabels := ispan(-6, 6, 2)*0.1
  end if 

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.02


  res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1,1/)
  res@xyMarkerSizes = (/0.03,0.03,0.03/)
  res@xyMarkerColors= (/"red","blue","black"/)


  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = (/"red","blue","black"/)
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  plot(k)  = gsn_csm_xy (wks, ispan(0,len-1,1), (/ssta_str_evo,ssta_weak_evo,ssta_nor_evo/), res)
 
  respl = True
  respl@gsLineColor = "gray"
  respl@gsLineDashPattern = 1.0

  x =  1.0*res@tmXBValues
  y = x

  if k.eq.0 then 
  y = 0.43*std_nino
  end if 

  if k.eq.1 then 
  y = -0.43*std_nino
  end if 
  
  if k.ne.2 then 
  plot@$stss(k)$ =  gsn_add_polyline(wks, plot(k),x, y, respl)
  end if



  ;;添加legend 
    gres = True
    if k.eq.0 then 
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    else
    gres@YPosPercent = 25.    ; expressed as %, 0->100, sets position of top border of legend 
    end if
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True

    textres@lgLabels = (/"Strong ISMR","Weak ISMR", "Normal ISMR"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
   

end do
pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,2/), pres)
end 