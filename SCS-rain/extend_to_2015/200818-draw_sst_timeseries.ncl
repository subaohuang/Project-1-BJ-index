
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201612

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1

  ;;
;;;region
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 170

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_annual = month_to_annual(ssta(:,{lat1:lat2},{lon1:lon2}), 1)
    ; ssta_annual =(/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)

    sst_mean = dim_avg_n_Wrap(ssta_annual(:,{0:15},{125:160}), (/1,2/))
;;;
   time = ispan(1979, 2016, 1)

;;;draw 
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "sst的时间序列-2016-philipian")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2016
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1980,2016,5)
  resxy@tmXBLabels = ispan(1980,2016,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2016,1)
  
  resxy@tmYLMode  =  "Explicit"
  resxy@tmYLValues = fspan(-0.4, 0.4,5)
  resxy@tmYLLabels = fspan(-0.4, 0.4, 5)
  resxy@trYMinF = -0.4
  resxy@trYMaxF = 0.4

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6




  resxy@xyLineThicknessF = 5
  resxy@xyLineColors       = "gray"

  plot = gsn_csm_xy(wks,time,sst_mean, resxy)


  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  lineres@lgLineColors = (/"gray"/) ; line colors
  ; lineres@lgLineThicknesses = (/ 3,3,3,3,3,3/)
  lineres@lgLineThicknesses = 5.
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
 

  textres = True
  textres@lgLabels = (/"SST"/)
  ; textres@lgLineLabelFontHeights = 0.005
  plot = simple_legend(wks,plot,gres,lineres,textres)   
draw(plot)
frame(wks)
end
