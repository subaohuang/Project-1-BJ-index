begin
;;;read cmip location trend and intension trend
;;;读入CMIP实验的结果
models1=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FGOALS-g2","FGOALS-s2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models1=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","FGOALS-s2","GFDL-CM3","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-MR","NorESM1-M"/)

 diri="/mnt/f/share/CMIP5/msf_history/"
 dim1=dimsizes(models1)
 print("dim="+dim1)
 ; width_trend_cmip=asciiread(diri+"his_msf_locationwidth_trend.ascii" ,dim1, "float")
 width_trend_cmip=asciiread(diri+"his_msf_locationxi_trend.ascii" ,dim1, "float")
 intension_trend_cmip=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim1, "float")
 
delete(diri)
;;;读入AMIP实验的结果
diri="/mnt/f/share/CMIP5/msf_AMIP/"
models2=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-s2","FGOALS-g2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)

  dim=dimsizes(models2)
  print("dim="+dim)
  ; width_trend_amip=asciiread(diri+"msf_locationwidth_trend.ascii" ,dim, "float")
  width_trend_amip=asciiread(diri+"msf_locationxi_trend.ascii" ,dim, "float")
  intension_trend_amip=asciiread(diri+"amip_msf_inten_trend.ascii" ,dim, "float")

  delete(dim)
  delete(diri)

;;;读入观测的结果
 rc_sst_rea_mean=0.34
 dim=6
 ; width_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_locationwidth_trend.ascii" ,dim, "float")
 width_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_locationxi_trend.ascii" ,dim, "float")
 intension_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_inten_trend.ascii" ,dim, "float")
;;;;



;;;平均的情况
width_trend_cmip_mean=avg(width_trend_cmip)
width_trend_amip_mean=avg(width_trend_amip)
width_trend_rea_mean=avg(width_trend_rea)


intension_trend_cmip_mean=avg(intension_trend_cmip)
intension_trend_amip_mean=avg(intension_trend_amip)
intension_trend_rea_mean=avg(intension_trend_rea)
;;;;



 ;;;画图，散点图
 ; rc1=regline(intension_trend_cmip,wid)
 ; line = rc*width_trend+rc@yintercept
 ; print("rc="+rc)

    ; mono = 1                            ; ascending=1 , descending=-1
    ; ii   =  dim_pqsort_n(width_trend,mono,0)
    ; xx   = width_trend(ii)                        ; ascending order 
    ; yy1   = sst_trend(ii)
printVarSummary(intension_trend_rea)
printVarSummary(width_trend_rea)



rc_cmip=regline(width_trend_cmip, intension_trend_cmip)
rc_amip=regline(width_trend_amip, intension_trend_amip)

line_cmip=rc_cmip*width_trend_cmip+rc_cmip@yintercept
line_amip=rc_amip*width_trend_amip+rc_amip@yintercept

r_cmip=escorc(width_trend_cmip, intension_trend_cmip)
r_amip=escorc(width_trend_amip, intension_trend_amip)
print("r_cmip="+r_cmip)
print("r_amip="+r_amip)
siglvl=0.99
  pr     = rtest(r_amip, dimsizes(width_trend_amip), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  pr     = rtest(r_cmip, dimsizes(width_trend_cmip), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)


r_rea=escorc(width_trend_rea, intension_trend_rea)

print("r_rea="+r_rea)
  pr     = rtest(r_rea, dimsizes(width_trend_rea), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/模式强度趋势和Wxi趋势的关系")
    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnMaximize=False
    res@xyMarkerColor =(/"white"/)

    res@xyMarkLineModes   = (/"Markers","lines"/)
    res@xyLineColor     = (/"black"/)          ; change line color
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    res@xyMarker         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.01      ; Marker size (default 0.01)

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF=0.8
    res@trYMinF=-0.8
    res@tmYLValues    = fspan(-0.8,0.8,9)
    res@tmYLLabels    = fspan(-0.8,0.8,9)
    
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

    res@vpXF=0.10
    res@vpYF=0.9
    res@vpHeightF=0.6
    res@vpWidthF=0.6
    
    ; res@tiMainString="Changes in the Walker Circulation"
    res@tiYAxisString="~F22~Trend in PWC intensity [10~S~9~N~Kg s~S~-1~N~/10yrs]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~Trend in ~F22~Width [~S~o~N~longitude/10yrs]"
    res@tiXAxisFontHeightF=0.012
    
    plot=gsn_csm_xy(wks,  width_trend_cmip, (/intension_trend_cmip,line_cmip/) ,res)
    res@xyMarkerColor =(/"white"/)
     res@xyLineColor     = (/"blue"/)          ; change line color
    
    plot1=gsn_csm_xy(wks,  width_trend_amip,(/intension_trend_amip,line_amip/),res)
    

    ;;;add 再分析资料的maker
    index_maker=(/2,15,5,6,7,8/)
    plot_rea=new(6,graphic)
    markers_rea=True
    markers_rea@gsMarkerSizeF=0.012
    markers_rea@gsMarkerThicknessF=2.0
    
    markers_rea@gsMarkerColor=(/"red"/)
    do i=0,5
      markers_rea@gsMarkerIndex=index_maker(i)
    plot_rea(i)=gsn_add_polymarker(wks,plot, width_trend_rea(i),intension_trend_rea(i), markers_rea)
    end do
    overlay(plot, plot1)
; 90~S~o~N~E

; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFontColor="black"
  txres@tfPolyDrawOrder="PostDraw"
  ; txres@txp
  labels1            = "  " + ispan(1,dim1,1)
  
  labels2           = (/"a","b","c","d","e","f","g","h","i","g","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)
  


  x1=conform(res@tmYLValues,0,0)
  y1=conform(res@tmXBValues,0,0)
  
  txid = gsn_add_text(wks,plot,labels1,width_trend_cmip,intension_trend_cmip,txres)
    txres@txFontColor="blue"
  printVarSummary(width_trend_amip)
  txid1 = gsn_add_text(wks,plot,labels2,width_trend_amip,intension_trend_amip,txres)

  
  polyline=True
  polyline@gsLineColor="black"
  polyline@gsLineThicknessF=2.0

  plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)

  plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

  Polymarker=True
  Polymarker@gsMarkerColor="black"
  Polymarker@gsMarkerSizeF=15.0
  Polymarker@gsMarkerIndex=16

  plotmark=gsn_add_polymarker(wks, plot, width_trend_cmip_mean, intension_trend_cmip_mean, Polymarker)
  Polymarker@gsMarkerColor="blue"
  plotmark2=gsn_add_polymarker(wks, plot, width_trend_amip_mean, intension_trend_amip_mean, Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark3=gsn_add_polymarker(wks, plot, width_trend_rea_mean, intension_trend_rea_mean, Polymarker)
 

;;;添加注释

labels=(/"1: ACCESS1-0","2: ACCESS1-3","3: BCC-CSM1-1","4: BCC-CSM1-1-M","5: CCSM4","6: CESM1-BGC","7: CMCC-CM","8: CMCC-CMS","9: CNRM-CM5",\
        "10: CSIRO-Mk3-6-0","11: FGOALS-g2","12: FGOALS-s2","13: FIO-ESM","14: GFDL-CM3","15: GFDL-ESM2G","16: GFDL-ESM2M",\
         "17: HadGEM2-AO","18: HadGEM2-CC","19: HadGEM2-ES","20: INMCM4","21: IPSL-CM5A-LR",\
    "22: IPSL-CM5A-MR","23: IPSL-CM5B-LR","24: MIROC5","25: MPI-ESM-LR","26: MPI-ESM-MR","27: MRI-CGCM3","28: NorESM1-M"/)

  xtxt=0.72
  ytxt = fspan(0.3,0.9,dim1)

  ytxt=ytxt(::-1)
   
  txres2               = True
  txres2@txFontHeightF = 0.010
  txres2@txFontColor="black"
  txres2@txJust="TopLeft"
   do i=0,dim1-1
      gsn_text_ndc (wks,labels(i),xtxt,ytxt(i),txres2)
    end do
    txres2@txFontThicknessF=3.0
    txres2@txFontHeightF = 0.015
    txres2@txFontColor="black"
     gsn_text_ndc (wks,"History",0.72,0.94,txres2)
     xtxt=0.90
      txres2@txFontColor="blue"
     gsn_text_ndc (wks,"AMIP",0.86,0.94,txres2)

labels2=(/"a: ACCESS1-0","b: ACCESS1-3","c: BCC-CSM1-1","d: BCC-CSM1-1-M",\
         "e: BNU-ESM","f: CanAM4","g: CCSM4","h: CMCC-CM","i: CNRM-CM5","j: CSIRO-Mk3-6-0","k: EC-EARTH","l: FGOALS-g2","m: FGOALS-s2",\
         "n: GFDL-CM3","o: GFDL-HIRAM-C180","p: HadGEM2-A","q: INMCM4","r: IPSL-CM5A-LR",\
    "s: IPSL-CM5A-MR","t: IPSL-CM5B-LR","u: MIROC5","v: MPI-ESM-LR","w: MPI-ESM-MR","x: MRI-AGCM3-2H","y: MRI-CGCM3","z: NorESM1-M"/)
  dim2=dimsizes(labels2)
  ytxt2=fspan(0.3, 0.9, 33)
  ytxt2=ytxt2(::-1)
  xtxt=0.86
  txres2@txFontHeightF = 0.010
  txres2@txFontColor="blue"
  do i=0,dim2-1
  gsn_text_ndc (wks,labels2(i),xtxt,ytxt2(i),txres2)
  end do

 ;;;;再分析的注释
  labels3=(/"20CR","ERA-Interium","JRA55","MERRA","NCEP1","NCEP2"/)
  xleg = 0.86
  xtxt_rea = 0.88  ; legend markers
  yleg = fspan(0.26, 0.38, 6)
  yleg=yleg(::-1)
  ytxt_rea = yleg  ; strings.

 do i=0,5
  gsres               = True
  gsres@gsMarkerIndex = index_maker(i)         ; Use filled dots for markers.
  gsres@gsMarkerColor      = "red"
  gsres@gsMarkerThicknessF = 0.015
  gsres@gsMarkerSizeF=10
   
; Add marker and text for the legend.
  gsn_polymarker_ndc(wks, xleg,yleg(i),gsres)
   txres2@txJust="CenterLeft"
   txres2@txFontColor="red"
  gsn_text_ndc      (wks,labels3(i),xtxt_rea,ytxt_rea(i),txres2)
end do

  txres2@txFontHeightF = 0.015
  txres2@txJust="TopLeft"
  gsn_text_ndc (wks,"Reanalysis",0.86,ytxt2(26),txres2)





;;;添加x轴下标
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
   gsn_text_ndc (wks,"shrinking",0.15,0.26,txres2)
   gsn_text_ndc (wks,"widening",0.65,0.26,txres2)
   
  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
  txres2@txJust="CenterCenter"
  txres2@txAngleF="90"
  plot_text=gsn_add_text(wks, plot, "strengthening",-7.5, 0.3, txres2)
  plot_text2=gsn_add_text(wks, plot, "weakening",-7.5, -0.3, txres2)
   

 ;;添加legend 
    gres = True
    gres@YPosPercent = 96.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 6    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    lineres@lgLineColors = (/"black","blue"/) ; line colors
    lineres@lgLineThicknesses = 2.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"r="+sprintf ("%4.2f",r_cmip),"r="+sprintf ("%4.2f",r_amip)/)

    plot = simple_legend(wks,plot,gres,lineres,textres)   

 draw(plot)
 frame(wks)

  end