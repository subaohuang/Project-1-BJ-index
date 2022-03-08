load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
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
function composite_mean(data1,data2,data3,data4,data5,data6,data7,data8,data9,yearind1)
begin 
    
    data_all = [/data1,data2,data3,data4,data5,data6,data7,data8,data9/]
    data_mean = new((/9,15/),float)

   do i=0,8
    data = data_all[i]
    data_mean(i,0:9) = dim_avg_n_Wrap(data(2:11,yearind1), 1)  ;;;M开始，MAM  JJASON  DJ
    data_mean(i,10:14) = dim_avg_n_Wrap(data(0:4,yearind1+1), 1)
   end do

return data_mean
end 


;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 187101
    ; endmon = 201912
    endmon = 201012
    startyear = 1871
    ; endyear = 2019
    endyear = 2010

;;读入HadiSST海温数据
    diri="/home/ys17-23/Extension/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019


    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;
    delete(date)
    delete(timeind)
    ;;
    
    ; diri = "/home/ys17-23/Extension/SODA_region/"
    ; fils1=systemfunc("ls "+diri+"soda_TEMP_1871-2010.nc")

    ; f1=addfile(fils1, "r")
    ; date:=cd_calendar(f1->TIME, 1)
    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    ; sst = f1->TEMP(timeind,0,:,:)

    ; ; sst=lonFlip(f1->sst(timeind,:,:))
    ; printVarSummary(sst)

    ; ssta=(/rmMonAnnCycTLL(sst)/)
    ; copy_VarCoords(sst, ssta)

    ; ;;detrend 
    ; ; ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ; ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&TIME),1),ssta,False,True,0)/)
    ; ;;;
    ; delete(date)
    ; delete(timeind)
;===========================================================================     
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_jja = nino34_djf
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

    nino34_djf = nino34_jja

    std_nino = dim_stddev(nino34_djf)

;;;


;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("../All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2010})  ;;;1871-2012

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2010-startyear+1,1),pre_jjasa,False,True,0)/)
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
    total_year = ispan(startyear+1,2009,1)   ;;;1872-2011
    DMI_ctrl_JJA_2 = DMI_ctrl_JJA(1:dimsizes(DMI_ctrl_JJA)-2)  ;;;1871-2011
    std = dim_stddev_Wrap(DMI_ctrl_JJA)


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


 nino34_djf_2 = nino34_djf(1:dimsizes(nino34_djf)-1-(endyear-2009))

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



    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)
  

;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


    case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

    case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

    case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]



    case = [/case1,case2,case3/]



    wksstring = (/"El Nino","La Nina","Normal"/)


    f4 = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/SODA/hb_soda.nc", "r")

    time := f4->time 
    date:=cd_calendar(time, 1)
    print(date(0))
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    dTdt_ct = f4->dTdt(timeind,:,:)

    UATXB_ct = f4->UATXB(timeind,:,:)  ;;;; ZA 256
    VATYB_ct = f4->VATYB(timeind,:,:)

    WATZB_ct = f4->WATZB(timeind,:,:)  ;;;EK

    UBTXA_ct = f4->UBTXA(timeind,:,:)  ;;MEAN zonal advection
    VBTYA_ct = f4->VBTYA(timeind,:,:)
    WBTZA_ct = f4->WBTZA(timeind,:,:) ;;; TH
    Q_ct = f4->Q(timeind,:,:)

    UATXA_ct = f4->UATXA(timeind,:,:)
    VATYA_ct = f4->VATYA(timeind,:,:)
    WATZA_ct = f4->WATZA(timeind,:,:)


    NLT = UATXA_ct
    NLT = (/UATXA_ct+VATYA_ct+WATZA_ct/)

;;;;
   
   ; lat1 = -10
   ; lat2 = 0
   ; lon1 = 270
   ; lon2 = 280

   ; lat1 = -5
   ; lat2 = 5
   ; lon1 = 210
   ; lon2 = 270

   ; lat1 = -5
   ; lat2 = 5
   ; lon1 = 150
   ; lon2 = 200


;     nino_type = (/"Nino1+2","Nino3","Nino3.4","Nino4"/)

; do j=0,3


;   if nino_type(j).eq."Nino1+2" then 
;    lat1 = -10
;    lat2 = 0
;    lon1 = 270
;    lon2 = 280
;    end if 

;   if nino_type(j).eq."Nino3" then 
;     lat1 = -5
;    lat2 = 5
;    lon1 = 210
;    lon2 = 270
;    end if 

;   if nino_type(j).eq."Nino3.4" then 
    lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
;    end if 

;   if nino_type(j).eq."Nino4" then 
;    lat1 = -5
;    lat2 = 5
;    lon1 = 150
;    lon2 = 200
;    end if 
  
  
;;;regional mean 
    dTdt_ct_nino3  =  dim_avg_n_Wrap(dTdt_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
 
    ; dTdt_ct_nino3  :=  dim_avg_n_Wrap(dsst(:,{lat1:lat2},{lon1:lon2}), (/1,2/))


    UATXB_ct_nino3  =  dim_avg_n_Wrap(UATXB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    VATYB_ct_nino3  =  dim_avg_n_Wrap(VATYB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    WATZB_ct_nino3 =  dim_avg_n_Wrap(VATYB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;;EK
    UBTXA_ct_nino3 = dim_avg_n_Wrap(UBTXA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    VBTYA_ct_nino3 = dim_avg_n_Wrap(VBTYA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    WBTZA_ct_nino3 = dim_avg_n_Wrap(WBTZA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    Q_ct_nino3 = dim_avg_n_Wrap(Q_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    NLT_nino3 = dim_avg_n_Wrap(NLT(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
;;;;

;;;split to 12 month 
    dTdt_ct_nino3_12 = split_12(dTdt_ct_nino3)    
    UATXB_ct_nino3_12 = split_12(UATXB_ct_nino3)
    VATYB_ct_nino3_12 = split_12(VATYB_ct_nino3)
    WATZB_ct_nino3_12 = split_12(WATZB_ct_nino3)

    UBTXA_ct_nino3_12 = split_12(UBTXA_ct_nino3)
    VBTYA_ct_nino3_12 = split_12(VBTYA_ct_nino3)
    WBTZA_ct_nino3_12 = split_12(WBTZA_ct_nino3)

    Q_ct_nino3_12 = split_12(Q_ct_nino3)
    NLT_nino3_12 = split_12(NLT_nino3)
;;;

  
  plot = new(3,graphic)
  midplot = plot
;;;;;

do k=0,2 

   case_need:= case[k]

   ind_ts1 := get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

   hb_mean := new((/3,9,15/), float)


   hb_mean_str = composite_mean(dTdt_ct_nino3_12,UATXB_ct_nino3_12,VATYB_ct_nino3_12,WATZB_ct_nino3_12,UBTXA_ct_nino3_12,\
                            VBTYA_ct_nino3_12,WBTZA_ct_nino3_12,Q_ct_nino3_12,NLT_nino3_12,ind_ts1)

   hb_mean_nor = composite_mean(dTdt_ct_nino3_12,UATXB_ct_nino3_12,VATYB_ct_nino3_12,WATZB_ct_nino3_12,UBTXA_ct_nino3_12,\
                            VBTYA_ct_nino3_12,WBTZA_ct_nino3_12,Q_ct_nino3_12,NLT_nino3_12,ind_tnor1)

   hb_mean_weak = composite_mean(dTdt_ct_nino3_12,UATXB_ct_nino3_12,VATYB_ct_nino3_12,WATZB_ct_nino3_12,UBTXA_ct_nino3_12,\
                            VBTYA_ct_nino3_12,WBTZA_ct_nino3_12,Q_ct_nino3_12,NLT_nino3_12,ind_tw1)

   hb_mean(0,:,:) =  hb_mean_str
   hb_mean(1,:,:) =  hb_mean_nor
   hb_mean(2,:,:) =  hb_mean_weak

   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,3:6), 2)   ;;;第三维0 是3月，6月是3
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,5:8), 2)   ;;;第三维0 是3月，6月是3
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,6:9), 2)   ;;SOND
   ; hb_mean_mamjja = dim_sum_n_Wrap(hb_mean(:,:,6:9), 2)   ;;SOND
   hb_mean_mamjja = dim_sum_n_Wrap(hb_mean(:,:,5:9), 2)   ;;ASON

   print(hb_mean_mamjja(0,0))
   print(sum(hb_mean_mamjja(0,1:)))
   ; exit()

;;**************************************************************
;画图
;;***************************************************************
     wks = gsn_open_wks("pdf",""+wksstring(k)+"-BAR-sum-JASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.3
     ; sres@vpXF = .15


     sres@trXMinF = 0.5
     sres@trXMaxF = 9.5
     sres@trYMinF = -1.2
     sres@trYMaxF = 1.2

     sres@gsnDraw = False
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.3           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = ispan(1, 9, 1)
     sres@tmXBLabels = (/"dTdt","ZA","VA","EK","MAU","MAV","TH","NHF","NDH"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "[~S~o~N~C~N~]" 
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
    

     str = ispan(1,2,1)
     ; color = (/"black","gray","yellow","green","pink","deeppink3","red","blue"/)
     color = (/"red"/)
     
    
     sres@gsnXYBarChartColors = color

     sres@gsnLeftStringFontHeightF = 0.02
     ; sres@gsnLeftString = "AMJJASOND mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino3 "
     sres@gsnLeftString = "JASON SUM Heat Budget for El Nino"
     ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(:,i),sres)   ; draw each time series

     plot1 = gsn_csm_xy(wks,ispan(1,9,1)-0.3,hb_mean_mamjja(0,:),sres)   ; draw each time series     
     color = (/"black"/)
     sres@gsnXYBarChartColors = color
     plot2 = gsn_csm_xy(wks,ispan(1,9,1),hb_mean_mamjja(1,:),sres)   ; draw each time series     

     sres@gsnXYBarChartColors = (/"blue"/)
     plot3 = gsn_csm_xy(wks,ispan(1,9,1)+0.3,hb_mean_mamjja(2,:),sres)   ; draw each time series     
    

     ; lbres                    = True          ; labelbar only resources
     ; lbres@vpWidthF           = 0.3           ; labelbar width
     ; lbres@vpHeightF          = 0.1           ; labelbar height
     ; lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
     ; lbres@lbFillColors       = (/"red","black"/)
     ; lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     ; lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
     ; lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     ; lbres@lbPerimOn          = False
     ; lbres@lgPerimColor    = "white"
     ; labels = (/"MSE","CpT"/)
     ; ; gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.33,lbres) ; draw right labelbar column
     ; gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.23,lbres) ; draw right labelbar column
  
     ; lbres@lbFillColors       := (/"green","blue"/)
     ; labels := (/"Lvq","Gz"/)
     ; ; gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.33,lbres) ; draw left labelbar column
     ; gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.23,lbres) ; draw left labelbar column
     
     ; ; pres = True
     ; ; gsn_panel(wks, plot, (/1,2/), pres)
   

     ; Markres = True
     ; Markres@gsMarkerColor = "gray"
     ; Markres@gsMarkerSizeF = 0.01
     ; Markres@gsMarkerIndex = 16
     ; plot@maker = gsn_add_polymarker(wks, plot(1),ispan(1,4,1), (/std_mse,std_A_mam,std_B_mam,std_C_mam/), Markres)
     ; draw(plot(1))
     ; frame(wks)

draw(plot1)
draw(plot2)
draw(plot3)
frame(wks)
end do 
; pres = True
; ; pres@gsnPanelLabelBar  = True
; ; pres@lbLabelStride=1
; gsn_panel(wks, plot, (/2,1/), pres)
end 