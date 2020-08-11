;;********************************************************
;;draw PC2 from MVEOF then compare it to 
;;timeseries of IPO index calculated from HadISST

;;lnshuheng 2020.01.09
;**********************************************************
begin
    ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    sst_annual = month_to_annual(ssta, 1)
    dim = dimsizes(sst_annual)
    
    ;;;;ipo index 
     sst_ipo = (/dim_avg_n_Wrap(sst_annual(:,{5:-5},{210:270}), (/1,2/))\
    -  (dim_avg_n_Wrap(sst_annual(:,{45:25},{150:210}), (/1,2/))+dim_avg_n_Wrap(sst_annual(:,{-25:-45},{150:210}), (/1,2/)))/2\
    /)
    sst_ipo_sm=runave_n_Wrap(sst_ipo, 3, 1, 0)
    ; sst_ipo_sm = sst_ipo
    sst_ipo_sm=(/dim_standardize_n_Wrap(sst_ipo_sm, 0, 0)/)
    
    sst_index = sst_ipo_sm

    sst_index_rc = regline(ispan(1, dimsizes(sst_index), 1), sst_index)

    sst_index = (/dtrend_msg_n(ispan(1,dimsizes(sst_index),1),sst_index,False,True,0)/)
    sst_index = (/-1*sst_index/)
  
    ;;;
    delete(timeind)
    delete(date)

    ;;;read MV- EOF PC1 data

    PC = asciiread("/home/ys17-19/lsh/Project/SCS-rain/sst_Q_pre_PC2_1958_2013.ascii", 56, "float")
    
    id1 = 1979-1958

    IPO = PC(id1:)
    printVarSummary(IPO)
    printVarSummary(sst_index)
    
    IPO=(/dim_standardize_n_Wrap(IPO, 0, 0)/)
    IPO = (/dtrend_msg_n(ispan(1,dimsizes(sst_index),1),IPO ,False,True,0)/)


  wks2   = gsn_open_wks("pdf", "IPO index和PC1进行比较-1979")
  plot2 = new(1,graphic)

  startyear = 1979
  endyear = 2013

  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False

  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = endyear
  resxy@trXMinF = startyear
  resxy@tmXBValues = ispan(startyear,endyear,5)
  resxy@tmXBLabels = ispan(startyear,endyear,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(startyear,endyear,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
  resxy@xyMonoLineThickness = True
  resxy@xyLineColors       =  (/"black","blue"/)          ; change line color
  resxy@xyMarkLineModes    =  (/"lines"/)
  resxy@xyMonoMarkLineMode = True


  time = ispan(startyear, endyear, 1)
  plot2 = gsn_csm_xy(wks2,time,(/IPO,sst_index/), resxy)
  draw(plot2)
  frame(wks2)
end 