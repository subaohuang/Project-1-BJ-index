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
function bw_bandpass(var)
 begin
    ca    = 9.0*12       ; band start (longer period)
    cb    = 4.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    opt   = False       ; use default options (time series of filtered
                        ; values will be returned)

    var = (/bw_bandpass_filter (var,fca,fcb,opt,0)/)  
 return var
 end 

;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    ; endmon = 201912
    endmon = 99912
    startyear = 1
    endyear = 999

;;;;
    diri_in = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"
    ; ;;;
    f4 = addfile(diri_in+"/hb/hb_picontrol.nc", "r")
    time := f4->time 
    date:=cd_calendar(time, 1)
    print(date(0))
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    dTdt = f4->dTdt(timeind,:,:)

    UATXB = f4->UATXB(timeind,:,:)  ;;;; ZA 256
    VATYB = f4->VATYB(timeind,:,:)
    WATZB = f4->WATZB(timeind,:,:)  ;;;EK

    UBTXA = f4->UBTXA(timeind,:,:)  ;;MEAN zonal advection
    VBTYA = f4->VBTYA(timeind,:,:)
    WBTZA = f4->WBTZA(timeind,:,:) ;;; TH
    Q= f4->Q(timeind,:,:)

    UATXA = f4->UATXA(timeind,:,:)
    VATYA = f4->VATYA(timeind,:,:)
    WATZA = f4->WATZA(timeind,:,:)

    NLT = UATXA
    NLT = (/UATXA+VATYA+WATZA/)
;;;;

   dim = dimsizes(dTdt)
;;;detrend 
   dTdt = (/dtrend_msg_n(ispan(1,dim(0),1),dTdt,False,True,0)/)
   UATXB = (/dtrend_msg_n(ispan(1,dim(0),1),UATXB,False,True,0)/)
   VATYB  = (/dtrend_msg_n(ispan(1,dim(0),1), VATYB ,False,True,0)/)
   WATZB  = (/dtrend_msg_n(ispan(1,dim(0),1), WATZB ,False,True,0)/)

   UBTXA = (/dtrend_msg_n(ispan(1,dim(0),1),UBTXA,False,True,0)/)
   VBTYA  = (/dtrend_msg_n(ispan(1,dim(0),1), VBTYA ,False,True,0)/)
   WBTZA  = (/dtrend_msg_n(ispan(1,dim(0),1),WBTZA ,False,True,0)/)
  
   Q = (/dtrend_msg_n(ispan(1,dim(0),1),Q,False,True,0)/)
   NLT = (/dtrend_msg_n(ispan(1,dim(0),1),NLT,False,True,0)/)

   dTdt = (/bw_bandpass(dTdt)/)

   UATXB = (/bw_bandpass(UATXB)/)
   VATYB = (/bw_bandpass(VATYB)/)
   WATZB= (/bw_bandpass(WATZB)/)

   UBTXA = (/bw_bandpass(UBTXA)/)
   VBTYA = (/bw_bandpass(VBTYA)/)
   WBTZA= (/bw_bandpass(WBTZA)/)

   Q = (/bw_bandpass(Q)/)
   
   NLT = (/bw_bandpass(NLT)/)

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile(diri_in+"piControl_monsoon_and_ensoyear_highpass_bw.nc", "r")

    nino_year = year_out->nino_year 
    lanina_year = year_out->lanina_year 
    normal_year = year_out->normal_year 

    monstryear = year_out->monstryear 
    monwekyear = year_out->monwekyear
    noryear = year_out->noryear
 
    wet_warm_year =  year_out->wet_warm_year
    wet_cold_year = year_out->wet_cold_year
    wet_norm_year = year_out->wet_norm_year

    dry_warm_year = year_out->dry_warm_year
    dry_cold_year = year_out->dry_cold_year
    dry_norm_year = year_out->dry_norm_year

    nor_warm_year = year_out->nor_warm_year
    nor_cold_year = year_out->nor_cold_year
    nor_norm_year = year_out->nor_norm_year
;;;;;


;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


    case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

    case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

    case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]
    case = [/case1,case2,case3/]
    wksstring = (/"El Nino","La Nina","Normal"/)

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
   ; lat1 = -10
   ; lat2 = 0
   ; lon1 = 270
   ; lon2 = 280
;    end if 

;   if nino_type(j).eq."Nino3" then 
   ;  lat1 = -5
   ; lat2 = 5
   ; lon1 = 210
   ; lon2 = 270
;    end if 

; ;   if nino_type(j).eq."Nino3.4" then 
   lat1 = -5
   lat2 = 5
   ; lon1 = 190
   ; lon2 = 240
   lon1 = 160
   lon2 = 270
; ;    end if 

;   if nino_type(j).eq."Nino4" then 
   ; lat1 = -5
   ; lat2 = 5
   ; lon1 = 150
   ; lon2 = 200
;    end if 
  
  
;;;regional mean 
    dTdt_nino3  =  dim_avg_n_Wrap(dTdt(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
 
    ; dTdt_nino3  :=  dim_avg_n_Wrap(dsst(:,{lat1:lat2},{lon1:lon2}), (/1,2/))


    UATXB_nino3  =  dim_avg_n_Wrap(UATXB(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    VATYB_nino3  =  dim_avg_n_Wrap(VATYB(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    WATZB_nino3 =  dim_avg_n_Wrap(WATZB(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;;EK
    UBTXA_nino3 = dim_avg_n_Wrap(UBTXA(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    VBTYA_nino3 = dim_avg_n_Wrap(VBTYA(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    WBTZA_nino3 = dim_avg_n_Wrap(WBTZA(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    Q_nino3 = dim_avg_n_Wrap(Q(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    NLT_nino3 = dim_avg_n_Wrap(NLT(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
;;;;

;;;split to 12 month 
    dTdt_nino3_12 = split_12(dTdt_nino3)    
    UATXB_nino3_12 = split_12(UATXB_nino3)
    VATYB_nino3_12 = split_12(VATYB_nino3)
    WATZB_nino3_12 = split_12(WATZB_nino3)

    UBTXA_nino3_12 = split_12(UBTXA_nino3)
    VBTYA_nino3_12 = split_12(VBTYA_nino3)
    WBTZA_nino3_12 = split_12(WBTZA_nino3)

    Q_nino3_12 = split_12(Q_nino3)
    NLT_nino3_12 = split_12(NLT_nino3)
;;;

  
  plot = new(3,graphic)
  midplot = plot

  ntilde  = "n~H-13V2F35~D~FV-2H3~"
;;;;;

do k=0,2 

   case_need:= case[k] ;;;  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

   ind_ts1 := get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

   hb_mean := new((/3,9,15/), float)


   hb_mean_str = composite_mean(dTdt_nino3_12,UATXB_nino3_12,VATYB_nino3_12,WATZB_nino3_12,UBTXA_nino3_12,\
                            VBTYA_nino3_12,WBTZA_nino3_12,Q_nino3_12,NLT_nino3_12,ind_ts1)

   hb_mean_nor = composite_mean(dTdt_nino3_12,UATXB_nino3_12,VATYB_nino3_12,WATZB_nino3_12,UBTXA_nino3_12,\
                            VBTYA_nino3_12,WBTZA_nino3_12,Q_nino3_12,NLT_nino3_12,ind_tnor1)

   hb_mean_weak = composite_mean(dTdt_nino3_12,UATXB_nino3_12,VATYB_nino3_12,WATZB_nino3_12,UBTXA_nino3_12,\
                            VBTYA_nino3_12,WBTZA_nino3_12,Q_nino3_12,NLT_nino3_12,ind_tw1)

   hb_mean(0,:,:) =  hb_mean_str
   hb_mean(1,:,:) =  hb_mean_nor
   hb_mean(2,:,:) =  hb_mean_weak

   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,3:6), 2)   ;;;第三维0 是3月，6月是3
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,5:8), 2)   ;;;第三维0 是3月，6月是3
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,6:9), 2)   ;;SOND
   ; hb_mean_mamjja = dim_sum_n_Wrap(hb_mean(:,:,6:9), 2)   ;;SOND
   ; hb_mean_mamjja = dim_sum_n_Wrap(hb_mean(:,:,5:9), 2)   ;;ASON 
                                                            ;;;MAMJJASONDJF
   
   hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,3:8), 2)   ;;JJASON
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,5:8), 2)   ;;JJASON  for ;;;;JASON
   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,4:7), 2)   ;;JJASON

   ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,0:9), 2)   ;;ASON
   

   print(hb_mean_mamjja(0,0))
   print(sum(hb_mean_mamjja(0,1:)))
   ; exit()

;;**************************************************************
;画图
;;***************************************************************
     wks = gsn_open_wks("pdf",""+wksstring(k)+"pIcontrol-BAR-sum-JASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf",""+wksstring(k)+"-BAR-sum-JASO-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf",""+wksstring(k)+"-BAR-sum-JJAS-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.3
     ; sres@vpXF = .15


     sres@trXMinF = 0.5
     sres@trXMaxF = 9.5

     ; sres@trYMinF = -1.2
     ; sres@trYMaxF = 1.2

     if k.ne.2 then 
     sres@trYMinF = -0.2
     sres@trYMaxF = 0.2
     else 
     sres@trYMinF = -0.1
     sres@trYMaxF = 0.1
     end if 

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
     sres@tiYAxisString = "[~S~o~N~C~N~/month]" 
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
     sres@gsnLeftString = "(a) Jun-Nov heat budget for El Ni"+ntilde+"o"
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