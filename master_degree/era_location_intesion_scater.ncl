;era数据计算的各强度指数和西边界的位置指数的相关性以及各个强度之间的相关性
function   xiangguan(intension,location)
local r1,prob,yes_no,siglvl
begin 
  r1=escorc_n(intension, location, 0, 0)
  prob=rtest(r1,dimsizes(intension),0)
  siglvl=0.05 
  yes_no = where(prob.lt.siglvl, 1.0, 0.0)
   ; print("500hpa以上msf_r1=" + r1 + "yes_no=" +yes_no)
   result=new(1,double)
   result(0)=r1
   ; result(1)=yes_no
return result
end 

function get_an_ao(sst)
begin
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    
    ssta_annual = month_to_annual(ssta_de, 1)
    ssta_annual!0 ="time"
    ssta_annual&time = ispan(1, 38, 1)
return ssta_annual
end 

 function  fitline(data1,data2)
 begin

    rc1=regCoef(data1,data2)
    line = rc1*data1+rc1@yintercept
  return line
  end
;************************************************************
begin 
 siglvl=0.05  ;检验值

;************************************************************
;;;read location dataset
;************************************************************
    ;;;read
    locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
    locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)
    locationxi_an = month_to_annual(locationxi, 1)
    ; locationxi_anu_std=dim_standardize(locationxi_anu, 1)
    locationxi_an!0="time"
    printVarSummary(locationxi_an)
    locationxi_an&time=ispan(1,38,1)
    ;;;提取数据
;;

;-----------------------------------------------------------------
;1.dslp 
;-----------------------------------------------------------------

   ;;read sur data
    fils = systemfunc("ls /mnt/f/share/era/data/data_series/sur_1.5/era_interim_sur1.5_"+"*"+".nc")
    f1 = addfiles(fils, "r")
    time = f1[:]->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
  
    u10 = short2flt(f1[:]->u10(timeind,{30:-30},:))
    slp = short2flt(f1[:]->msl(timeind,{30:-30},:))
    slp = (/slp/100/) ;;change to hpa
    v10 = short2flt(f1[:]->v10(timeind,{30:-30},:))
    
    delete(fils)
    delete(f1)
    delete(date)
    delete(timeind)
    delete(time)
    ;;;
    

    ; slp_a = get_an_ao(slp)
    ; u10_a = get_an_ao(u10)

    du10 = dim_avg_n_Wrap(u10(:,{-5:5},{150:210}), (/1,2/))
    
    dslp = dim_avg_n_Wrap(slp(:,{-5:5},{200:280}),(/1,2/)) - dim_avg_n_Wrap(slp(:,{-5:5},{80:160}),(/1,2/))
    copy_VarCoords(slp(:,1,1), dslp)
    
     
    du10 = (/dtrend_msg_n(ispan(1,456,1),du10,False,True,0)/)
    dslp = (/dtrend_msg_n(ispan(1,456,1),dslp,False,True,0)/)
 

    du10_an = month_to_annual(du10, 1)
    dslp_an = month_to_annual(dslp, 1)

    du10_an = (/-1*du10_an/)

    




    ;******************************************************************************************
    ;STRF
    ;**********************************************************************************************
  
    fils = systemfunc("ls /mnt/f/share/rea_msf/msf_time_era_global.nc")
    f1 = addfile(fils, "r")
    time = f1->time
    date = cd_calendar(time, 1)
    timeind = ind(date.le.201612.and.date.ge.197901)
    msf_r = f1->msf_time({100:1000},:,timeind)
    msf_r = (/msf_r/10^9/)
    
    msf_r = (/dtrend_msg_n(ispan(1,456,1),msf_r,False,True,2)/)

    STRF = dim_avg_n_Wrap(msf_r(:,{150:210},:),(/0,1/))
    copy_VarCoords(msf_r(1,1,:), STRF)
    ; STRF_a = get_an_ao(STRF)
    STRF_an = month_to_annual(STRF, 1)
    delete(date)
    delete(timeind)
    delete(time)
    ;******************************************************************************************
    ;OMEGA
    ;**********************************************************************************************
    fils6 = systemfunc("ls /mnt/f/share/era/data/data_series/w_1.5/era_interim_moda_"+"*"+".nc")
    f6 = addfiles(fils6, "r")
    date=cd_calendar(f6[:]->time, 1)
    varname = getfilevarnames(f6[1])
    timeind=ind(date.le.201612.and.date.ge.197901)
    w_r =short2flt(f6[:]->w(timeind,{500},{5:-5},:))
    
    w_r =  (/dtrend_msg_n(ispan(1,456,1),w_r,False,True,0)/)
    w_d = dim_avg_n_Wrap(w_r(:,:,{200:240}),(/1,2/)) - dim_avg_n_Wrap(w_r(:,:,{120:160}),(/1,2/))
    copy_VarCoords(w_r(:,1,1), w_d)
    ; w_d_a = get_an_ao(w_d)
    printVarSummary(w_d)
    
    w_d_an = month_to_annual(w_d, 1)
    w_d_an = (/10*w_d_an/)
    ;;;;

 ;*****************************************************************
 ;各个指数的相关性  intension_msf,slp,u10_mean,intension_ud,STRF
 ;**************************************************************
   result=new(4,double)
   result(0)=xiangguan(locationxi_an,dslp_an)
   result(1)=xiangguan(locationxi_an,du10_an)
   result(2)=xiangguan(locationxi_an,w_d_an)
   result(3)=xiangguan(locationxi_an,STRF_an)



   do i=0,3
    print("result"+i+"=" + result(i) )
   end do
  ;*********************************************************************
  ;画强度与位置的散点图，表明情况
  ;************************************************************************
 

    line1=fitline(locationxi_an,dslp_an)
    line2=fitline(locationxi_an,du10_an)
    line3=fitline(locationxi_an,w_d_an)
    line4=fitline(locationxi_an,STRF_an)
    

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/scatter_intension&location")

    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@xyLineThicknesses = (/ 1.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"Blue"  ,"black"/)          ; change line color
    resxy@xyMarkerColor = "Blue"
    resxy@xyMarkLineModes   = (/"Markers","lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    resxy@xyMarkers         =  16                      ; choose type of marker  
    resxy@xyMarkerSizeF     = 0.005        ; Marker size (default 0.01)
     
    resxy@trXMinF = 140
    resxy@trXMaxF = 182

    resxy@gsnLeftString="(a)dSLP"
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01

    resxy@vpHeightF=0.3
    resxy@vpWidthF=0.6
    ; resxy@vpXF=0.15
    ; resxy@vpYF=0.95

    resxy@tmXBLabelFontAspectF=1.0
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用

    resxy@tmXBMinorLengthF=0.003
    resxy@tmXBMinorOutwardLengthF=0.003

    resxy@tmYLMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    resxy@tmYLLabelFont=30
    resxy@tmYLLabelFontHeightF=0.008
    resxy@tmYLMinorOn=False
    resxy@tmXTOn=False
    resxy@tmYROn=False


    plot=new(4,graphic)


    plot(0)=gsn_csm_xy(wks,locationxi_an,(/dslp_an,line1/),resxy)
    resxy@gsnLeftString="(b)-U10"
    ; resxy@gsnRightString=""+u10_data(3,1)+""

    plot(1)=gsn_csm_xy(wks,locationxi_an,(/du10_an,line2/),resxy)
    resxy@gsnLeftString="(c)omega"

    plot(2)=gsn_csm_xy(wks,locationxi_an,(/w_d_an,line3/),resxy)
    resxy@gsnLeftString="(d)STRF"
    ; resxy@tmYLMode="Explicit"    
    ; resxy@tmYLValues=fspan(0.4, 1.0, 7)
    ; resxy@tmYLLabels=fspan(0.4,1.0,7)
    plot(3)=gsn_csm_xy(wks,locationxi_an,(/STRF_an,line4/),resxy)


   

    txid_tr = new(4, graphic)
    amid_tr = new(4, graphic)
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.01

    amres_tr                  = True
    amres_tr@amParallelPosF   =  0.45   ; This is the right edge of the plot.
    amres_tr@amOrthogonalPosF =  -0.35    ; This is the top edge of the plot.

    r=result
    letter = sprintf("%5.3f",r)
do i = 0, 3
;---Create text strings
    tr_label = "("  +"r="+ letter(i) + ")" 

    txres@txBackgroundFillColor = "transparent"
    txres@txPerimOn = False
    txid_tr(i) = gsn_create_text(wks, tr_label, txres)

;---Attach text strings to plot
    amid_tr(i) = gsn_add_annotation(plot(i), txid_tr(i), amres_tr)

end do

  pres = True
  pres@gsnMaximize=True
  gsn_panel(wks, plot, (/2,2/), pres)
    ; maximize_output(wks,pres)  ; calls draw and frame for you 
    
 end
