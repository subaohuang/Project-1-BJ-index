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
 ;;;time span
    startmon =198001
    endmon = 202012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;


;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
    delete(date)
    delete(timeind)

    sst_ct = ssta
    sst_ct=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)


    sst_ct_djf = month_to_season(sst_ct, "DJF")
    sst_ct_djf(0:dimsizes(sst_ct_djf&time)-2,:,:) = sst_ct_djf(1:dimsizes(sst_ct_djf&time)-1,:,:)
    sst_ct_djf(dimsizes(sst_ct_djf&time)-1,:,:) = sst_ct(dimsizes(sst_ct&time)-1,:,:)

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
    enino_year_ct_or = sst_year(ind(nino3_4_ct.ge.0.5*std_ct_nino3_4))
    ; ; enino_year_ct = enino_year_ct_or(1:dimsizes(enino_year_ct_or)-2)
    
    ; ; enino_year_wnp = enino_year_wnp_or(1:dimsizes(enino_year_wnp_or)-2)

    print(enino_year_ct_or)
    ; exit()
    ; print(nino3_4_ct)

 
    ; enino_year_1 =(/1983,1987,1988,1992,1995,1998/) -1
    ; enino_year_2 =(/2003,2007,2010,2015,2016/) -1

    ;;;for lanina
    enino_year_1 =(/1989,1999/)-1
    enino_year_2 =(/2000,2008,2011,2012/) -1
    ; enino_year_1 = enino_year_ct_or(0:5)

    ; enino_year_2 = enino_year_ct_or(6:)
  ; enino_year_2 := (/2015,2015/)

    enino_year_ct = get1Dindex(sst_year,enino_year_1)
    enino_year_wnp = get1Dindex(sst_year,enino_year_2)




    ; ely=[1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016];
    enino_year = [/enino_year_ct,enino_year_wnp/]

    hb_mean = new((/2,9,15/),float)

    f4 = addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/B2000/heat_bud/obs/hb_godas.nc", "r")

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


    nino_type = (/"Nino1+2","Nino3","Nino3.4","Nino4"/)

do j=0,3


  if nino_type(j).eq."Nino1+2" then 
   lat1 = -10
   lat2 = 0
   lon1 = 270
   lon2 = 280
   end if 

  if nino_type(j).eq."Nino3" then 
    lat1 = -5
   lat2 = 5
   lon1 = 210
   lon2 = 270
   end if 

  if nino_type(j).eq."Nino3.4" then 
    lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
   end if 

  if nino_type(j).eq."Nino4" then 
   lat1 = -5
   lat2 = 5
   lon1 = 150
   lon2 = 200
   end if 
  
  
;;;regional mean 
    dTdt_ct_nino3  =  dim_avg_n_Wrap(dTdt_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
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


;;;composite mean 
do i=0,1  
    hb_mean_ct = composite_mean(dTdt_ct_nino3_12,UATXB_ct_nino3_12,VATYB_ct_nino3_12,WATZB_ct_nino3_12,UBTXA_ct_nino3_12,\
                            VBTYA_ct_nino3_12,WBTZA_ct_nino3_12,Q_ct_nino3_12,NLT_nino3_12,enino_year[i])

    hb_mean(i,:,:) = hb_mean_ct 
end do
  ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,0:5), 2)
  ; hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,1:9), 2)
  hb_mean_mamjja = dim_avg_n_Wrap(hb_mean(:,:,3:8), 2)


  print("**********************************")
  
  print("DT="+hb_mean_mamjja(0,0))
  print("SUM="+sum(hb_mean_mamjja(0,1:8)))
  

  plot = new(2,graphic)
  midplot = plot
   

;;**************************************************************
;画图
;;***************************************************************

     ; wks = gsn_open_wks("pdf","OBS-NINO3-M2D-BAR-MAMJJA-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf","NOAA-YEAR-NINO1+2-M2D-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf","NOAA-YEAR-NINO3-M2D-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf","LAnina-NOAA-YEAR-"+nino_type(j)+"-M2D-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     ; wks = gsn_open_wks("pdf","El Nino-NOAA-YEAR-"+nino_type(j)+"-M2D-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file
     wks = gsn_open_wks("pdf","La Nina-NOAA-YEAR-"+nino_type(j)+"-M2D-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15


     sres@trXMinF = 0.5
     sres@trXMaxF = 9.5
     sres@trYMinF = -0.5
     sres@trYMaxF = 0.5

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
     color = (/"black"/)
     
    
     sres@gsnXYBarChartColors = color

     sres@gsnLeftStringFontHeightF = 0.02
     ; sres@gsnLeftString = "AMJJASOND mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino3 "
     sres@gsnLeftString = "GODAS JJASON mean Heat Budget"+" "+nino_type(j)
     ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(:,i),sres)   ; draw each time series

     plot1 = gsn_csm_xy(wks,ispan(1,9,1)-0.15,hb_mean_mamjja(0,:),sres)   ; draw each time series     
     color = (/"red"/)
     sres@gsnXYBarChartColors = color
     plot2 = gsn_csm_xy(wks,ispan(1,9,1)+0.15,hb_mean_mamjja(1,:),sres)   ; draw each time series     
    

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
frame(wks)
end do 
; pres = True
; ; pres@gsnPanelLabelBar  = True
; ; pres@lbLabelStride=1
; gsn_panel(wks, plot, (/2,1/), pres)
end 