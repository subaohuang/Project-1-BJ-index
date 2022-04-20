begin
;;;read cmip location trend and sst trend
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","FGOALS-s2","GFDL-CM3","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-MR","NorESM1-M"/)
 diri="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 print("dim="+dim)
 width_trend=asciiread(diri+"his_msf_locationwidth_trend.ascii" ,dim, "float")
 sst_trend=asciiread("/mnt/f/share/CMIP5/sst_cha/his_sstca_trend.ascii" ,dim, "float")
 printVarSummary(width_trend)
 print(width_trend)
 print(sst_trend)
 rc_sst_rea_mean=0.32
 width_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_locationwidth_trend.ascii" ,dim, "float")

 rc_width_trend_rea_mean=avg(width_trend_rea)
 print("rc_width_trend_rea_mean="+rc_width_trend_rea_mean)

 ;;;画图，散点图
 rc=regline(width_trend,sst_trend)
 line = rc*width_trend+rc@yintercept
 print("rc="+rc)
 r=escorc(width_trend, sst_trend)
 print("r="+r)

    mono = 1                            ; ascending=1 , descending=-1
    ii   =  dim_pqsort_n(width_trend,mono,0)
    xx   = width_trend(ii)                        ; ascending order 
    yy1   = sst_trend(ii)


width_trend_mean=avg(width_trend)
sst_trend_mean=avg(sst_trend)

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/cmip模式海温差趋势和宽度趋势的关系")
    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@xyLineColors     = (/"red"/)          ; change line color
    res@xyMarkerColor ="white"
    res@xyMarkLineModes   = (/"Markers","Lines"/)
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.01      ; Marker size (default 0.01)

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=0.4
    res@trYMinF=-0.4
    res@tmYLValues    = fspan(-0.4,0.4,9)
    res@tmYLLabels    = fspan(-0.4,0.4,9)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=8
    res@trXMinF=-8
    res@tmXBValues    = ispan(-8,8,2)
    res@tmXBLabels    = ispan(-8,8,2)

    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.008  ;调整刻度的长度，注意2个命令要同时使用
    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8

    res@tmYLLabelFontHeightF=0.015
    res@vpXF=0.15
    res@vpYF=0.9
    res@vpHeightF=0.6
    res@vpWidthF=0.6
    
    res@tiMainString="Changes in the Walker Circulation"
    res@tiYAxisString="~F22~Trend in ~F5~D~F22~SST [k/10yrs]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~Trend in ~F22~Width [~S~o~N~longitude/10yrs]"
    res@tiXAxisFontHeightF=0.012
    
   plot=gsn_csm_xy(wks, width_trend, (/sst_trend,line/), res)
; 90~S~o~N~E

; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFontColor="blue"
  txres@tfPolyDrawOrder="PostDraw"
  ; txres@txp
  labels            = "  " + ispan(1,dim,1)
  
  x1=conform(res@tmYLValues,0,0)
  y1=conform(res@tmXBValues,0,0)
  
  txid = gsn_add_text(wks,plot,labels,width_trend,sst_trend,txres)
  
  polyline=True
  polyline@gsLineColor="black"
  polyline@gsLineThicknessF=2.0

  plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)

  plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

  Polymarker=True
  Polymarker@gsMarkerColor="blue"
  Polymarker@gsMarkerSizeF=15.0
  Polymarker@gsMarkerIndex=16

  plotmark=gsn_add_polymarker(wks, plot, width_trend_mean, sst_trend_mean, Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark2=gsn_add_polymarker(wks, plot, rc_width_trend_rea_mean, rc_sst_rea_mean, Polymarker)

  xtxt=0.8
  ytxt = fspan(0.3,0.9,dim)



labels=(/"1: ACCESS1-0","2: ACCESS1-3","3: BCC-CSM1-1","4: BCC-CSM1-1-M","5: CCSM4","6: CESM1-BGC","7: CMCC-CM","8: CMCC-CMS","9: CNRM-CM5",\
        "10: CSIRO-Mk3-6-0","11: FIO-ESM","12: GFDL-CM3","13: GFDL-ESM2G","14: GFDL-ESM2M",\
         "15: HadGEM2-AO","16: HadGEM2-CC","17: HadGEM2-ES","18: INMCM4","19: IPSL-CM5A-LR",\
    "20: IPSL-CM5A-MR","21: IPSL-CM5B-LR","22: MIROC5","23: MPI-ESM-LR","24: MPI-ESM-MR","25: MRI-CGCM3","26: NorESM1-M"/)


  txres2               = True
  txres2@txFontHeightF = 0.015
  txres2@txFontColor="blue"
  txres2@txJust="TopLeft"
   do i=0,dim-1
      gsn_text_ndc (wks,labels(i),xtxt,ytxt(dim-1-i),txres2)
    end do
  txres2@txFontThicknessF=3.0
    gsn_text_ndc (wks,"CMIP5",xtxt,0.95,txres2)
;;;添加x轴下标
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
   gsn_text_ndc (wks,"shrinking",0.2,0.26,txres2)
   gsn_text_ndc (wks,"widening",0.7,0.26,txres2)
   


  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
  txres2@txAngleF="90"
  plot_text=gsn_add_text(wks, plot, "strengthening",-7, 0.3, txres2)
  plot_text2=gsn_add_text(wks, plot, "weakening",-7, -0.3, txres2)
   
  
  ;;;
  
  txres2@txFontColor="black"
  txres2@txAngleF="0"
  plot_text3=gsn_add_text(wks, plot, "Correlation="+sprintf ("%4.2f",r),-4, 0.37, txres2)



   ; ;;legend 
   ;  gres = True
   ;  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   ;                      ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
   ;  gres@XPosPercent = 70    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

   ;  lineres = True
   ;  lineres@lgLineColors = (/"black"/) ; line colors
   ;  lineres@lgLineThicknesses = 3.0                       ; line thicknesses
   ;  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

   ;  textres = True
   ;  textres@lgLabels = "Correlation="+sprintf ("%4.2f",r)

   ;  plot = simple_legend(wks,plot,gres,lineres,textres)   

 draw(plot)
 frame(wks)

  end