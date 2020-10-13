load "/mnt/d/ncl/Dot-Matrix-master/WCHAR.ncl" ; All related function tar into this lib 
begin
;;;read cmip location trend and sst trend
;;;read cmip location trend and sst trend
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CanESM2","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; ;FGOALS-g2和s2没有海温数据l;;同时添加了models=(/"CanESM2",GISS-E2-R","MIROC-ESM-CHEM","NorESM1-ME"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CCSM4","CanESM2","CESM1-BGC","CESM1-CAM5","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

 models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M",\
    "BNU-ESM","CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
    "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
    "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME"/)

    ;;;
    diri="/mnt/f/share/CMIP5/msf_history/"
    dim=dimsizes(models)

    width_trend=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim, "float")
    sst_trend=asciiread("/mnt/f/share/CMIP5/sst_cha/his_sstca_trend.ascii" ,dim, "float")
    printVarSummary(width_trend)
    print(width_trend)
    print(sst_trend)
    rc_sst_rea_mean=0.22
    rc_width_trend_rea_mean=0.31
    print("rc_width_trend_rea_mean="+rc_width_trend_rea_mean)

    ;;;画图，散点图,两条拟合线
    
    rc=regline(width_trend,sst_trend)

    rc1=regline(width_trend(0:17),sst_trend(0:17))
    rc2=regline(width_trend(18:),sst_trend(18:))


    width_trend_new=new(37,typeof(width_trend))
    width_trend_new(1:35)=width_trend
    width_trend_new(0)=-0.2
    width_trend_new(36)=0.2

    line= rc*width_trend_new+rc@yintercept
    ; line1 = rc1*width_trend(0:17)+rc1@yintercept
    ; line2= rc2*width_trend(18:)+rc2@yintercept
    r=escorc(width_trend, sst_trend)
    R2=r^2
    ; r1=escorc(width_trend(0:17), sst_trend(0:17))
    ; r2=escorc(width_trend(18:), sst_trend(18:))
    ; print("r1="+r1)
    ; print("r2="+r2)
    ; R21=r1^2
    ; R22=r2^2


 
    width_trend_mean=avg(width_trend)
    sst_trend_mean=avg(sst_trend)

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/scs_fig_意见图_多模式.6")
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
    
    res@tmXTOn=True
    res@tmYROn=True
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXTMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=0.4
    res@trXMinF=-0.4
    res@tmXBValues    = fspan(-0.4,0.4,9)
    res@tmXBLabels    = fspan(-0.4,0.4,9)
     res@tmXTValues    = fspan(-0.4,0.4,9)
    res@tmXTLabels    = fspan(-0.4,0.4,9)

    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用
    
    res@tmXTMajorLengthF=0.008        ;調整刻度的长度
    res@tmXTMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用
    
    res@tmYLMajorLengthF=0.008    
    res@tmYLMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用    ;調整刻度的长度

    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8
    res@tmBorderThicknessF=2.0

    res@tmYLLabelFontHeightF=0.015
    res@vpXF=0.13
    res@vpYF=0.9
    res@vpHeightF=0.6
    res@vpWidthF=0.58
    
    ; res@tiMainString="Changes in the Walker Circulation"

    res@tiYAxisString="~F5~D~F22~SST        /(kdecade~S~-1~N~)"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~PWCI         /(10~S~9~N~Kg S~S~-1~N~decade~S~-1~N~)"
    res@tiXAxisFontHeightF=0.012
    
   plot=gsn_csm_xy(wks, width_trend, (/sst_trend/), res)
  
    res_poly=True
    res_poly@gsLineThicknessF=3.0
    res_poly@gsLineColor="black"

    plot1=gsn_add_polyline(wks, plot, width_trend_new, line, res_poly)
    ; plot1=gsn_add_polyline(wks, plot, width_trend(0:17), line1, res_poly)
    ; res_poly@gsLineColor="blue"
    ; plot2=gsn_add_polyline(wks, plot, width_trend(18:), line2, res_poly)


; 90~S~o~N~E

; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.02
  txres@txFontColor="black"
  txres@tfPolyDrawOrder="PostDraw"
  ; txres@txp
  ; labels            = "  " + ispan(1,dim,1)
  ; labels=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","M"/)
  labels= array_append_record ((/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/),""+ispan(1, 17, 1)+"",0)
  printVarSummary(labels)
                               
  x1=conform(res@tmYLValues,0,0)
  y1=conform(res@tmXBValues,0,0)
  
  txid = gsn_add_text(wks,plot,labels(0:17),width_trend(0:17),sst_trend(0:17),txres)
    txres@txFontColor="blue"
  txid1 = gsn_add_text(wks,plot,labels(18:),width_trend(18:),sst_trend(18:),txres)

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

  xtxt=0.74
  ytxt = fspan(0.3,0.9,18)
  ytxt=ytxt(::-1)



; labels=(/"1: ACCESS1-0","2: ACCESS1-3","3: BCC-CSM1-1","4: BCC-CSM1-1-M","5: CCSM4","6: CESM1-BGC","7: CMCC-CM","8: CMCC-CMS","9: CNRM-CM5",\
;         "10: CSIRO-Mk3-6-0","13: FIO-ESM","14: GFDL-CM3","15: GFDL-ESM2G","16: GFDL-ESM2M",\
;          "17: HadGEM2-AO","18: HadGEM2-CC","19: HadGEM2-ES","20: INMCM4","21: IPSL-CM5A-LR",\
;     "22: IPSL-CM5A-MR","23: IPSL-CM5B-LR","24: MIROC5","25: MPI-ESM-LR","26: MPI-ESM-MR","27: MRI-CGCM3","28: NorESM1-M"/)

; alfpe= (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D"/)
; alfpe=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)

alfpe= array_append_record ((/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/),""+ispan(1, 17, 1)+"",0)


labels = alfpe+ ":" +models


  txres2               = True
  txres2@txFontHeightF = 0.01
  txres2@txFontColor="black"
  txres2@txJust="TopLeft"

    do i=0,dim-18
     gsn_text_ndc (wks,labels(i),xtxt,ytxt(i),txres2)
    end do
     xtxt=0.86
  txres2@txFontColor="blue"
   do i=dim-17,dim-1
     gsn_text_ndc (wks,labels(i),xtxt,ytxt(i-18),txres2)
    end do

   txres2@txFontColor="black"
   txres2@txFontThicknessF=3.0
    gsn_text_ndc (wks,"CMIP5",0.825,0.95,txres2)
   

     ;;;;添加中文注释
    fontres=True
    fontres@FontQuality="High"
    fontres@fontAngleF=90
    fontres@fontSize=10
    fontres@fontNDC=True

    fontX(wks,plot, "趋势", 0.077,0.578, fontres)  
    fontres@fontAngleF=0
    fontX(wks,plot, "趋势", 0.40,0.255, fontres)
  

    ; drawNDCGrid(wks)
  
  ;;;
                                    
  ; txres2@txFontColor="black"
  ; txres2@txAngleF="0"
  ; plot_text3=gsn_add_text(wks, plot, "Correlation="+sprintf ("%4.2f",r),-4, 0.37, txres2)



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
    textres@lgLabels = (/"R~S~2~N~="+sprintf("%4.2f",R2)/)
    plot = simple_legend(wks,plot,gres,lineres,textres)


 draw(plot)
 frame(wks)

  end