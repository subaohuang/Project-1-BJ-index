begin
;;;read cmip location trend and sst trend
;;;read cmip location trend and sst trend
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CanESM2","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)


 diri="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 print("dim="+dim)
 width_trend=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim, "float")
 sst_trend=asciiread("/mnt/f/share/CMIP5/sst_cha/his_sstca_trend.ascii" ,dim, "float")
 printVarSummary(width_trend)
 print(width_trend)
 print(sst_trend)
 rc_sst_rea_mean=0.22
 rc_width_trend_rea_mean=0.31
 print("rc_width_trend_rea_mean="+rc_width_trend_rea_mean)

 ;;;画图，散点图
 rc=regline(width_trend,sst_trend)


 width_trend_new=new(20,typeof(width_trend))
 width_trend_new(1:18)=width_trend
 width_trend_new(0)=-0.2
 width_trend_new(19)=0.2
  line = rc*width_trend_new+rc@yintercept
  
 print("rc="+rc)
 r=escorc(width_trend, sst_trend)
 print("r="+r)
 R2=r^2

    mono = 1                            ; ascending=1 , descending=-1
    ii   =  dim_pqsort_n(width_trend,mono,0)
    xx   = width_trend(ii)                        ; ascending order 
    yy1   = sst_trend(ii)
     rc2=regline(xx,yy1)

 
width_trend_mean=avg(width_trend)
sst_trend_mean=avg(sst_trend)

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/cmip模式海温差趋势和强度趋势的关系_2008_有AMIP的模式")
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
    res@trYMaxF=0.3
    res@trYMinF=-0.3
    res@tmYLValues    = fspan(-0.3,0.3,7)
    res@tmYLLabels    = fspan(-0.3,0.3,7)
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=0.4
    res@trXMinF=-0.4
    res@tmXBValues    = fspan(-0.4,0.4,9)
    res@tmXBLabels    = fspan(-0.4,0.4,9)

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
    res@tmBorderThicknessF=4.0

    res@tmYLLabelFontHeightF=0.015
    res@vpXF=0.15
    res@vpYF=0.9
    res@vpHeightF=0.6
    res@vpWidthF=0.6
    
    res@tiMainString="Changes in the Walker Circulation"
    res@tiYAxisString="~F22~Trend in ~F5~D~F22~SST [k/10yrs]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~Trend in ~F22~PWCI [10~S~9~N~Kg S~S~-1~N~/10yrs]   "

    res@tiXAxisFontHeightF=0.012
    
   plot=gsn_csm_xy(wks, width_trend, (/sst_trend/), res)
  
res_poly=True
res_poly@gsLineThicknessF=3.0
res_poly@gsLineColor="black"
plot1=gsn_add_polyline(wks, plot, width_trend_new, line, res_poly)

; 90~S~o~N~E

; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.02
  txres@txFontColor="black"
  txres@tfPolyDrawOrder="PostDraw"
  ; txres@txp
  ; labels            = "  " + ispan(1,dim,1)
  ; labels=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","M"/)
  labels=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)
  x1=conform(res@tmYLValues,0,0)
  y1=conform(res@tmXBValues,0,0)
  
  txid = gsn_add_text(wks,plot,labels,width_trend,sst_trend,txres)
  
  polyline=True
  polyline@gsLineColor="black"
  polyline@gsLineThicknessF=2.0

  plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)

  plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

  Polymarker=True
  Polymarker@gsMarkerColor="black"
  Polymarker@gsMarkerSizeF=15.0
  Polymarker@gsMarkerIndex=16

  plotmark=gsn_add_polymarker(wks, plot, width_trend_mean, sst_trend_mean, Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark2=gsn_add_polymarker(wks, plot, rc_width_trend_rea_mean, rc_sst_rea_mean, Polymarker)

  xtxt=0.8
  ytxt = fspan(0.3,0.9,dim)




; labels=(/"1: ACCESS1-0","2: ACCESS1-3","3: BCC-CSM1-1","4: BCC-CSM1-1-M","5: CCSM4","6: CESM1-BGC","7: CMCC-CM","8: CMCC-CMS","9: CNRM-CM5",\
;         "10: CSIRO-Mk3-6-0","13: FIO-ESM","14: GFDL-CM3","15: GFDL-ESM2G","16: GFDL-ESM2M",\
;          "17: HadGEM2-AO","18: HadGEM2-CC","19: HadGEM2-ES","20: INMCM4","21: IPSL-CM5A-LR",\
;     "22: IPSL-CM5A-MR","23: IPSL-CM5B-LR","24: MIROC5","25: MPI-ESM-LR","26: MPI-ESM-MR","27: MRI-CGCM3","28: NorESM1-M"/)

; alfpe= (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D"/)
alfpe=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)
labels = alfpe+ ":" +models

  txres2               = True
  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
  txres2@txJust="TopLeft"
   do i=0,dim-1
     gsn_text_ndc (wks,labels(i),xtxt,ytxt(dim-1-i),txres2)
    end do
  txres2@txFontThicknessF=3.0
    gsn_text_ndc (wks,"CMIP5",xtxt,0.95,txres2)
;;;添加x轴下标
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
   gsn_text_ndc (wks,"weakening",0.2,0.26,txres2)
   gsn_text_ndc (wks,"strengthening",0.7,0.26,txres2)
   


  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
  txres2@txAngleF="90"
   gsn_text_ndc (wks,"weakening",0.08,0.40,txres2)
   gsn_text_ndc (wks,"strengthening",0.08,0.80,txres2)
   
  
  ;;;
                                    
  txres2@txFontColor="black"
  txres2@txAngleF="0"
  plot_text3=gsn_add_text(wks, plot, "Correlation="+sprintf ("%4.2f",r),-4, 0.37, txres2)



   ;;legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    lineres@lgLineColors = (/"black"/) ; line colors
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = "regression"+"("+"R~S~2~N~="+sprintf("%4.2f",R2)+")"

    plot = simple_legend(wks,plot,gres,lineres,textres)   

 draw(plot)
 frame(wks)

  end