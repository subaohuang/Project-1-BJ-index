begin

;;;读入AMIP实验的结果
diri="/mnt/f/share/CMIP5/msf_AMIP/"
models2=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-s2","FGOALS-g2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)

  dim=dimsizes(models2)
  print("dim="+dim)
  width_trend_amip=asciiread(diri+"msf_locationwidth_trend.ascii" ,dim, "float")
  intension_trend_amip=asciiread(diri+"amip_msf_inten_trend.ascii" ,dim, "float")
  
  width_trend_alt_amip=asciiread(diri+"msf_locationwidth_alt_trend.ascii" ,dim, "float")
  intension_trend_alt_amip=asciiread(diri+"amip_msf_inten_alt_trend.ascii" ,dim, "float")

  width_trend_id_amip=asciiread(diri+"msf_locationxi_trend.ascii" ,dim, "float")
  intension_trend_id_amip=asciiread(diri+"amip_msf_inten_id_trend.ascii" ,dim, "float")



  delete(dim)
  delete(diri)

;;;读入观测的结果
 ; rc_sst_rea_mean=0.34
 dim=6
 width_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_locationwidth_trend.ascii" ,dim, "float")
 intension_trend_rea=asciiread("/mnt/f/share/rea_msf/rea_inten_trend.ascii" ,dim, "float")

 width_trend_alt_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_alt_locationwidth_trend.ascii" ,dim, "float")
 intension_trend_alt_rea=asciiread("/mnt/f/share/rea_msf/rea_inten_alt_trend.ascii" ,dim, "float")

 width_trend_id_rea=asciiread("/mnt/f/share/rea_msf/rea_msf_locationxi_trend.ascii" ,dim, "float")
 intension_trend_id_rea=asciiread("/mnt/f/share/rea_msf/rea_inten_id_trend.ascii" ,dim, "float")

;;;;



;;;平均的情况
width_trend_amip_mean=avg(width_trend_amip)
width_trend_rea_mean=avg(width_trend_rea)


intension_trend_amip_mean=avg(intension_trend_amip)
intension_trend_rea_mean=avg(intension_trend_rea)
;;;
width_trend_alt_amip_mean=avg(width_trend_alt_amip)
width_trend_alt_rea_mean=avg(width_trend_alt_rea)

intension_trend_alt_amip_mean=avg(intension_trend_alt_amip)
intension_trend_alt_rea_mean=avg(intension_trend_alt_rea)
;;;
width_trend_id_amip_mean=avg(width_trend_id_amip)
width_trend_id_rea_mean=avg(width_trend_id_rea)

intension_trend_id_amip_mean=avg(intension_trend_id_amip)
intension_trend_id_rea_mean=avg(intension_trend_id_rea)
;;;;

printVarSummary(intension_trend_rea)
printVarSummary(width_trend_rea)

width_rea_amip=new((/32/),typeof(width_trend_amip))
width_rea_alt_amip=new((/32/),typeof(width_trend_amip))
width_rea_id_amip=new((/32/),typeof(width_trend_amip))

width_rea_amip(0:25)=width_trend_amip
width_rea_amip(26:)=width_trend_rea
width_rea_alt_amip(0:25)=width_trend_alt_amip 
width_rea_alt_amip(26:)=width_trend_alt_rea
width_rea_id_amip(0:25)=width_trend_id_amip
width_rea_id_amip(26:)=width_trend_id_rea

intension_rea_amip=new((/32/),typeof(intension_trend_amip))
intension_rea_alt_amip=new((/32/),typeof(intension_trend_amip))
intension_rea_id_amip=new((/32/),typeof(intension_trend_amip))

intension_rea_amip(0:25)=intension_trend_amip
intension_rea_amip(26:)=intension_trend_rea
intension_rea_alt_amip(0:25)=intension_trend_alt_amip 
intension_rea_alt_amip(26:)=intension_trend_alt_rea
intension_rea_id_amip(0:25)=intension_trend_id_amip
intension_rea_id_amip(26:)=intension_trend_id_rea





rc_id=regline(width_rea_id_amip, intension_rea_id_amip)
rc_alt=regline(width_rea_alt_amip, intension_rea_alt_amip)
rc_pa=regline(width_rea_amip, intension_rea_amip)


line_id=rc_id*width_rea_id_amip+rc_id@yintercept
line_alt=rc_alt*width_rea_alt_amip+rc_alt@yintercept
line_pa=rc_pa*width_rea_amip+rc_pa@yintercept


r_id=escorc(width_rea_id_amip, intension_rea_id_amip)
r_alt=escorc(width_rea_alt_amip, intension_rea_alt_amip)
r_pa=escorc(width_rea_amip, intension_rea_amip)




print("r_id="+r_id)
print("r_pa="+r_pa)
print("r_alt="+r_alt)

siglvl=0.99
  pr     = rtest(r_id, dimsizes(width_rea_id_amip), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  pr     = rtest(r_pa, dimsizes(width_rea_amip), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)


  pr     = rtest(r_alt, dimsizes(width_rea_alt_amip), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)




  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/三个环流模式强度趋势和Width趋势的关系")
    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnMaximize=False
    res@xyMarkerColor =(/"white"/)
    res@xyMarkLineMode  = (/"Markers"/)
    
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

    ; res@vpXF=0.10
    ; res@vpYF=0.9
    ; res@vpHeightF=0.6
    ; res@vpWidthF=0.6
    
    res@tiYAxisString="~F22~Trend in PWC intensity [10~S~9~N~Kg s~S~-1~N~/10yrs]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~Trend in ~F22~Width [~S~o~N~longitude/10yrs]"
    res@tiXAxisFontHeightF=0.012
    ;;;

    plot=new(3,graphic)
    plot1=plot


    plot(0)=gsn_csm_xy(wks,  width_trend_id_amip ,intension_trend_id_amip,res)
    res@xyLineColor     = (/"blue"/)          ; change line color
    res@xyMarkLineModes   = (/"lines"/)
    plot1(0)=gsn_csm_xy(wks,  width_rea_id_amip, (/line_id/) ,res)
   ;;id
    delete(res@xyMarkLineModes)
     res@xyMarkLineMode   = (/"Markers"/)
    plot(1)=gsn_csm_xy(wks,  width_trend_amip ,intension_trend_amip,res)
    res@xyLineColor     = (/"blue"/)          ; change line color
    res@xyMarkLineModes   = (/"lines"/)
    plot1(1)=gsn_csm_xy(wks,  width_rea_amip, (/line_pa/) ,res)
    ;;alt
      res@xyMarkLineMode   = (/"Markers"/)
    plot(2)=gsn_csm_xy(wks,  width_trend_alt_amip ,intension_trend_alt_amip,res)
    res@xyMarkLineMode   = (/"lines"/)
    res@xyLineColor     = (/"blue"/)          ; change line color
    plot1(2)=gsn_csm_xy(wks,  width_rea_alt_amip, (/line_alt/) ,res)

  

    do i=0,2
  overlay(plot(i), plot1(i))
end do
 str=unique_string("dum")
 ;;;add 再分析资料的maker
    index_maker=(/2,15,5,6,7,8/)
    plot_rea=new(1,graphic)
    markers_rea=True
    markers_rea@gsMarkerSizeF=0.012
    markers_rea@gsMarkerThicknessF=2.0
    markers_rea@gsMarkerColor=(/"red"/)




  j=0
  do i=0,5
      markers_rea@gsMarkerIndex=index_maker(i)
    plot_rea@$unique_string("dum")$=gsn_add_polymarker(wks,plot(j), width_trend_id_rea(i),intension_trend_id_rea(i), markers_rea)
    end do
    j=1
 do i=0,5
      markers_rea@gsMarkerIndex=index_maker(i)
    plot_rea@$unique_string("dum")$=gsn_add_polymarker(wks,plot(j), width_trend_rea(i),intension_trend_rea(i), markers_rea)
    end do
    j=2
 do i=0,5
      markers_rea@gsMarkerIndex=index_maker(i)
    plot_rea@$unique_string("dum")$=gsn_add_polymarker(wks,plot(j), width_trend_alt_rea(i),intension_trend_alt_rea(i), markers_rea)
    end do



; Attach some text strings next to the markers.
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFontColor="black"
  txres@tfPolyDrawOrder="PostDraw"
  ; txres@txp

  
  labels2           = (/"a","b","c","d","e","f","g","h","i","g","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)
  


  x1=conform(res@tmYLValues,0,0)
  y1=conform(res@tmXBValues,0,0)
  ; txid1=new(3,graphic)
  
  txid1 = gsn_add_text(wks,plot(0),labels2,width_trend_id_amip,intension_trend_id_amip,txres)
  txid2 = gsn_add_text(wks,plot(1),labels2,width_trend_amip,intension_trend_amip,txres)
  txid3 = gsn_add_text(wks,plot(2),labels2,width_trend_alt_amip,intension_trend_alt_amip,txres)

  
  ; polyline=True
  ; polyline@gsLineColor="black"
  ; polyline@gsLineThicknessF=2.0

  ; plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)
  ; plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

  Polymarker=True
  Polymarker@gsMarkerColor="black"
  Polymarker@gsMarkerSizeF=15.0
  Polymarker@gsMarkerIndex=16

  plotmark2=gsn_add_polymarker(wks, plot(0), width_trend_id_amip_mean, intension_trend_id_amip_mean , Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark3=gsn_add_polymarker(wks, plot(0), width_trend_id_rea_mean, intension_trend_id_rea_mean, Polymarker)

  plotmark4=gsn_add_polymarker(wks, plot(1), width_trend_amip_mean, intension_trend_amip_mean, Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark5=gsn_add_polymarker(wks, plot(1), width_trend_rea_mean, intension_trend_rea_mean, Polymarker)
  
   plotmark6=gsn_add_polymarker(wks, plot(2), width_trend_alt_amip_mean, intension_trend_alt_amip_mean, Polymarker)
  Polymarker@gsMarkerColor="red"
  plotmark7=gsn_add_polymarker(wks, plot(2), width_trend_alt_rea_mean, intension_trend_alt_rea_mean, Polymarker)

;;;添加注释

   
  txres2               = True
  txres2@txFontHeightF = 0.010
  txres2@txFontColor="black"
  txres2@txJust="TopLeft"

  txres2@txFontThicknessF=3.0
  txres2@txFontHeightF = 0.015
  txres2@txFontColor="black"
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
  plot_text=new(1,graphic)

  do i=0,2
  plot_text@$unique_string("str1")$=gsn_add_text(wks, plot(i), "strengthening",-7.5, 0.3, txres2)
  plot_text@$unique_string("str2")$ = gsn_add_text(wks, plot(i), "weakening",-7.5, -0.3, txres2)
  end do

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
    r_amip=new(3,typeof(r_id))
    r_amip(0)=r_id
    r_amip(1)=r_pa
    r_amip(2)=r_alt
do i=0,2
  textres@lgLabels = (/sprintf ("%4.2f",r_amip(i))/)
    plot(i) = simple_legend(wks,plot(i),gres,lineres,textres)   
end do 

   pres=True

   gsn_panel(wks, plot, (/2,2/), pres)

  end