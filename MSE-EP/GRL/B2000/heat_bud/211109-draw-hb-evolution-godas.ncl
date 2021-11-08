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
function composite_mean(data1,data2,data3,data4,data5,data6,data7,data8,yearind1)
begin 
    
    data_all = [/data1,data2,data3,data4,data5,data6,data7,data8/]
    data_mean = new((/8,14/),float)

   do i=0,7 
    data = data_all[i]
    data_mean(i,0:9) = dim_avg_n_Wrap(data(2:11,yearind1), 1)
    data_mean(i,10:13) = dim_avg_n_Wrap(data(0:3,yearind1+1), 1)
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

    ; enino_year_ct_or = ind(nino3_4_ct.ge.std_ct_nino3)
    ; ; enino_year_ct = enino_year_ct_or(1:dimsizes(enino_year_ct_or)-2)
    
    ; ; enino_year_wnp = enino_year_wnp_or(1:dimsizes(enino_year_wnp_or)-2)

    ; print(enino_year_ct)
    ; print(nino3_4_ct)


    enino_year_1 =(/1980,1983,1987,1988,1992,1995,1998/)

    enino_year_2 =(/2003,2005,2007,2010,2015,2016/)

   sst_year = ispan(1980,2020,1)

   enino_year_ct = get1Dindex(sst_year,enino_year_1)
   enino_year_wnp = get1Dindex(sst_year,enino_year_2)


    ; ely=[1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016];


   enino_year = [/enino_year_ct,enino_year_wnp/]
 


  hb_mean = new((/2,8,14/),float)


do i=0,1
   
    f4 = addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/B2000/heat_bud/obs/hb_godas.nc", "r")


    time := f4->time 
    time = (/time-1/)
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
;;;;
   
   lat1 = -10
   lat2 = 0
   lon1 = 270
   lon2 = 280
  
;;;regional mean 
    dTdt_ct_nino3  =  dim_avg_n_Wrap(dTdt_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    UATXB_ct_nino3  =  dim_avg_n_Wrap(UATXB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    VATYB_ct_nino3  =  dim_avg_n_Wrap(VATYB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    WATZB_ct_nino3 =  dim_avg_n_Wrap(VATYB_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;;EK
    UBTXA_ct_nino3 = dim_avg_n_Wrap(UBTXA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    VBTYA_ct_nino3 = dim_avg_n_Wrap(VBTYA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    WBTZA_ct_nino3 = dim_avg_n_Wrap(WBTZA_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
    Q_ct_nino3 = dim_avg_n_Wrap(Q_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/)) ;;MEAN zonal advection
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
;;;


;;;composite mean 
   
    hb_mean_ct = composite_mean(dTdt_ct_nino3_12,UATXB_ct_nino3_12,VATYB_ct_nino3_12,WATZB_ct_nino3_12,UBTXA_ct_nino3_12,\
                                VBTYA_ct_nino3_12,WBTZA_ct_nino3_12,Q_ct_nino3_12,enino_year[i])
;;;

    hb_mean(i,:,:) = hb_mean_ct 
end do


  plot = new(2,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "NINO1+2_Heat-Budget-evolution-for-OBS-two-state")
   

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
    res@trXMaxF        = 13



    res@gsnYRefLine = 0.0


    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@gsnLeftStringFontHeightF   = 0.02
    res@gsnRightStringFontHeightF  = 0.02

    res@tmXBMode = "Explicit"
    res@tmXBValues = ispan(0, 13, 1)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
    res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A"/)
    res@tmYLMode = "Explicit"

    ; if k.ne.0 then 
    res@trYMinF        = -1.0
    res@trYMaxF        = 1.0

    res@tmYLValues := fspan(-1.0, 1.0, 11)
    res@tmYLLabels := fspan(-1.0, 1.0, 11)
    ; else
    ;  res@trYMinF        = -1.2
    ; res@trYMaxF        = 1.2

    ; res@tmYLValues := fspan(-1.2, 1.2, 13)
    ; res@tmYLLabels := fspan(-1.2, 1.2, 13)
    ; end if 

    res@tiYAxisString = "[~S~o~N~C~N~/month]" 
    res@tiYAxisFontHeightF=0.02

    res@tiXAxisString = "Month"
    res@tiXAxisFontHeightF = 0.02


    res@xyMarkLineMode = (/"MarkLines"/)
    res@xyMonoMarkLineMode = True
    res@xyMonoMarker = True
    res@xyMarker = (/1/)
    res@xyMonoMarkerSize = True
    res@xyMarkerSizeF = (/0.03/)

    ; res@xyMarkerColors= (/"firebrick","dimgray","deepskyblue3","deeppink3","darkslategray4",\
    ;                           "darkseagreen3","darkorchid3","dodgerblue4"/)


     res@xyMarkerColors = (/"black","gray","yellow","green","pink","deeppink3","red","blue"/)

    res@xyMonoLineThickness = True
    res@xyLineThicknessF = 2.0
    

    res@xyLineColors       = res@xyMarkerColors 
    res@xyMonoDashPattern = True
    res@xyDashPattern     = (/1/)
    res@lgPerimOn = False
    res@lgLabelFontHeightF = 0.018

    res@gsnLeftString = "(a) Heat budget terms for CTRL"

    
    plot(0) = gsn_csm_xy(wks, ispan(0,13,1), hb_mean(0,:,:), res)

     res@gsnLeftString = "(b) Heat budget terms for ATL warming"

    plot(1) = gsn_csm_xy(wks, ispan(0,13,1), hb_mean(1,:,:), res)
   
  ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
 ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgMonoDashIndex = True
    lineres@lgDashIndex = res@xyDashPattern
    textres = True

    textres@lgLabels = (/"DTdt","ZAU","ZAV","EK","MAU","MAV","TH","NHF"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   do i=0,1
    midplot(i) = simple_legend(wks,plot(i),gres,lineres,textres)   
   end do

pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,1/), pres)
end 