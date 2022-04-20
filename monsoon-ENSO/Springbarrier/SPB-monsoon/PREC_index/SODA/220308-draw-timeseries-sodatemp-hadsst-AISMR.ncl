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


;===========================================================   
    diri = "/home/ys17-23/Extension/SODA_region/"
    fils1=systemfunc("ls "+diri+"soda_TEMP_1871-2010.nc")

    f1:=addfile(fils1, "r")
    date:=cd_calendar(f1->TIME, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst_soda = f1->TEMP(timeind,0,:,:)

    ssta_soda = (/rmMonAnnCycTLL(sst_soda)/)
    copy_VarCoords(sst_soda, ssta_soda)

    ;;detrend 
    ssta_soda = (/dtrend_msg_n(ispan(1,dimsizes(ssta_soda&TIME),1),ssta_soda,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
;=======================================================================



    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")
    ; nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    ; nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 
    nino34_12_r =  split_12(nino34)
    ; nino34_jja = nino34_djf
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)
    ; nino34_djf = nino34_jja

    ; std_nino = dim_stddev(nino34_djf)
    

    nino34_soda =  dim_avg_n_Wrap(ssta_soda(:,{5:-5},{190:240}), (/1,2/))

    nino34_soda_djf = month_to_season(nino34_soda, "DJF")
    ; nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    ; nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 
    nino34_12_r_soda =  split_12(nino34_soda)
    ; nino34_jja = nino34_djf
    nino34_jja_soda = (/dim_avg_n_Wrap(nino34_12_r_soda(5:8,:), 0)/)
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
    DMI_ctrl_JJA = (/dim_standardize(DMI_ctrl_JJA, 1)/)
;;;
   
   nino34_jja = (/runave(nino34_jja, 3, 1)/)
   nino34_jja_soda = (/runave(nino34_jja_soda, 3, 1)/)
   DMI_ctrl_JJA = (/runave(DMI_ctrl_JJA, 3, 1)/)



x =  ispan(1871,2010,1)
plot = new(3,graphic)
wks = gsn_open_wks("PDF", "HadISST-Sodatemp-AISMR-index-timeseries")
;;======================================================================
;;draw
;;================================================================


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


  ; res@trXMinF        = 0
  ; res@trXMaxF        = 190

  res@gsnYRefLine = 0.0

  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  ; res@tmXBMode = "Explicit"
  ; res@tmXBValues = ispan(0, 17, 1)
  ; ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
   res@tmYLMode = "Explicit"

   res@trYMinF        = -2.0
   res@trYMaxF        = 2.0

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Years"
  res@tiXAxisFontHeightF = 0.02

  ; res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  ; res@xyMarkers = (/1,1,1/)
  ; res@xyMarkerSizes = (/0.03,0.03,0.03/)
  ; res@xyMarkerColors= (/"red","blue","black"/)


  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = (/"red","blue","black"/)
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  k =0
  plot(k)  = gsn_csm_xy (wks, x, (/nino34_jja,nino34_jja_soda,DMI_ctrl_JJA/), res)


  ; respl = True
  ; respl@gsLineColor = "gray"
  ; respl@gsLineDashPattern = 1.0

  ; x =  1.0*res@tmXBValues
  ; y = x

  ; if k.eq.0 then 
  ; y = 0.43*std_nino
  ; end if 

  ; if k.eq.1 then 
  ; y = -0.43*std_nino
  ; end if 
  
  ; if k.ne.2 then 
  ; plot@$stss(k)$ =  gsn_add_polyline(wks, plot(k),x, y, respl)
  ; end if


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
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True

    textres@lgLabels = (/"HADISST","SODATEMP", "AISMR"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
    midplot = plot
    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
  
  draw(plot(0))
  frame(wks)
; pres = True
; ; pres@gsnPanelLabelBar  = True
; ; pres@lbLabelStride=1
; gsn_panel(wks, plot, (/2,2/), pres)
end 
