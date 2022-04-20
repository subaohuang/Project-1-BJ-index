
;;************************************
;annual cycle of WL and SST
;**************************************

function annual_cycle(u10_inten)
local u10_inten_month, nmos,nmo
begin
    u10_inten_month = new(12,typeof(u10_inten))
    
    nmos = 12
    ntim = 456
    do nmo=0,nmos-1
    u10_inten_month(nmo) = dim_avg( u10_inten(nmo:ntim-1:nmos) )
    end do
return u10_inten_month
end

begin 
;;---------------------------------------
;;read data
;----------------------------------------
    locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
    locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)
    locationxi!0 = "time"
    locationxi&time = ispan(1,456,1)

    locationxi_month = new(12,typeof(locationxi))
    nmos = 12
    ntim = 456

   do nmo=0,nmos-1
      locationxi_month(nmo) = dim_avg( locationxi(nmo:ntim-1:nmos) )
   end do
   
;;;
;;;read sst and calculate correlation
    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    time = f2->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    sst=lonFlip(f2->sst(timeind,{30:-30},:))

    sst = (/dtrend_msg_n(ispan(1,456,1),sst,False,True,0)/)
    
    sst_clm   = clmMonTLL(sst)
    nino3 = dim_avg_n_Wrap( sst_clm(:,{-5:5},{210:270}) , (/1,2/))
 
    ; nino3  =  dim_avg_n_Wrap(sst_clm(:,{-5:5},{80:160}), (/1,2/))-dim_avg_n_Wrap(sst_clm(:,{-5:5},{100:280}), (/1,2/)) 
    delete([/time,date,timeind/])
;;;

;;---------------------------------------
;;read  wind data
;----------------------------------------
    fils = systemfunc("ls /mnt/f/share/era/data/data_series/sur_1.5/era_interim_sur1.5_"+"*"+".nc")
    f3 = addfiles(fils, "r")

    time=f3[:]->time
    date=cd_calendar(time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    u10 = short2flt(f3[:]->u10(tind,:,:))
    slp = short2flt(f3[:]->msl(tind,:,:))
    printVarSummary(u10)


    u10_inten = dim_avg_n_Wrap(u10(:,{-5:5},{150:240}),(/1,2/))
    u10_inten = (/-1*u10_inten/)
    u10_inten_month = annual_cycle(u10_inten)

    slp_inten = dim_avg_n_Wrap(slp(:,{-5:5},{200:280}),(/1,2/)) - \
                  dim_avg_n_Wrap(slp(:,{-5:5},{80:160}),(/1,2/))
    copy_VarCoords(slp(:,1,1), slp_inten)
                
    slp_inten_month = annual_cycle(slp_inten)



;; corr 
   
; ;*************************************************** 超前滞后相关
;     mxlag    = 12
;     x_Lead_y = esccr(nino3_com(2,:),u10_inten_com(2,:), mxlag)
;     y_Lead_x = esccr(locationxi_com(2,:),nino3_com(2,:), mxlag)    ; switch the order of the series

;     ccr = new ( 2*mxlag+1, double)    
;     ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
;     ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
;     x1 = ispan(-12,12,1)     
;     wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/不同ENSO相位下位置和强度相关性")
;     res=True
;     res@gsnXRefLine = 0
;     plot1=gsn_csm_xy(wks1,x1,ccr,res)
; ;;;

x = ispan(1,12,1)

;;---------------------------------------
;;draw_piture
;----------------------------------------

;;;
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/slp和NINO3区的年循环")

    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@xyLineThicknesses= (/ 2.0/)          ; make second line thicker
    resxy@xyLineColor      = (/"red"/)          ; change line color
    resxy@xyMarkerColor = "red"
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode   = (/"MarkLines"/) 
    

    resxy@xyMarkers         =  (/16/)                    ; choose type of marker  
    resxy@xyMarkerSizeF     = 0.01        ; Marker size (default 0.01)
    
  
    ; resxy@gsnRightString=""+slp_data(3,1)+""
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01

    resxy@trXMaxF = 12
    resxy@trXMinF = 1
    
    resxy@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    resxy@tmXBValues=ispan(1,12,1)
    resxy@tmXBLabels = ""+x

    
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
 
    

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = .2                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.3                ; move units down
    resxy@pmLegendWidthF         = 0.15                ; Change width and
    resxy@pmLegendHeightF        = 0.18                ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .01                ; label font height
    resxy@xyExplicitLegendLabels = (/"nino3"/)         ; create explicit labels


    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.6
    ; resxy@vpXF=0.15
    ; resxy@vpYF=0.95
    res2= resxy
    res2@xyLineColor      = (/"black"/)          ; change line color
    res2@xyMarkerColor = "black"
    res2@xyExplicitLegendLabels = (/"WL"/)         ; create explicit labels
    
    ; plot = gsn_csm_xy2(wks, ispan(1,12, 1), nino3, locationxi_month, resxy, res2)
     plot = gsn_csm_xy2(wks, ispan(1,12, 1), nino3, slp_inten_month, resxy, res2)

    draw(plot)
    frame(wks)
   end