function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 
function split_12_4d(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2),dim(3)/),float)
    do i=0,11   
    data_12(i,:,:,:,:) = data(i::12,:,:,:)
    end do 
return data_12
end 
begin
  yearstart = 198001
  yearend = 202012
  


  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360



  f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
  MSE = f1->MSE

; ;;;;read data
    var= (/"omega.mon.mean.nc","shum.mon.mean.nc","rhum.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)     
    w = f1->omega(timeind,{1000:300},:,:)
    
    w = (/-50*w/)

    q = f2->shum(timeind,{1000:300},:,:)

    q= q/1000.0

    rh = f3->rhum(timeind,{1000:300},:,:)
;;;


;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    sst_r=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst_r)

    sst_o = sst_r(:,::-1,:)

    delete(date)
    delete(timeind)
;;;

;;;;
    yo = fspan(-90.0, 90.0, 73)
    xo = fspan(0.0, 357.5, 144)
    sst = linint2_Wrap(sst_o&longitude, sst_o&latitude, sst_o, True, xo, yo, 0)

    delete(sst&LATITUDE)
    delete(sst&LONGITUDE)
    sst!1 = "lat"
    sst!2 = "lon"
    sst&lat = yo
    sst&lon = xo
;;;;

;;;;
  ssta = sst
  ssta=(/rmMonAnnCycTLL(sst)/)
  ssta_de = ssta
  ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta_de&time),1),ssta,False,True,0)/)
;;;;



;;;vertical mean  
  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)
   
  level = 850
  q_col =  dim_avg_wgt_n_Wrap(q, dz, 1, 1)
   ; q_col = q(:,{level},:,:)
  
  q_col = (/rmMonAnnCycTLL(q_col)/)
  q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
  
  w_col =  dim_avg_wgt_n_Wrap(w, dz, 1, 1)
    ; w_col = w(:,{level},:,:)
  w_col = (/rmMonAnnCycTLL(w_col)/)
  w_col = (/dtrend_msg_n(ispan(1,dimsizes(w_col&time),1),w_col,False,True,0)/)

  rh_col =  dim_avg_wgt_n_Wrap(rh, dz, 1, 1)
  ; rh_col = rh(:,{level},:,:)
  rh_col = (/rmMonAnnCycTLL(rh_col)/)
  rh_col = (/dtrend_msg_n(ispan(1,dimsizes(rh_col&time),1),rh_col,False,True,0)/)
;;;;

  yearnum = 2020-1980 +1
;;;


;;;;annual mean 
    q_12 = split_12_4d(q)
    q_annual = q(0:yearnum-1,:,:,:)
    q_annual = (/dim_avg_n_Wrap(q_12(0:4,:,:,:,:), 0)/)

    rh_12 = split_12_4d(rh)
    rh_annual = rh(0:yearnum-1,:,:,:)
    rh_annual = (/dim_avg_n_Wrap(rh_12(0:4,:,:,:,:), 0)/)

    w_12 = split_12_4d(w)
    w_annual = w(0:yearnum-1,:,:,:)
    w_annual = (/dim_avg_n_Wrap(w_12(0:4,:,:,:,:), 0)/)


    sst_12 = split_12(sst)
    sst_annual = sst(0:yearnum-1,:,:)
    sst_annual = (/dim_avg_n_Wrap(sst_12(0:4,:,:,:), 0)/)

    ; corr = escorc_n(q_col, p,0, 0)
    corr = escorc_n(rh_annual, w_annual,0, 0)
    ; corr = escorc_n(q_col_annual, w_col_annual,0, 0)
    copy_VarCoords(q_annual(1,:,:,:), corr)

    corr_mean = dim_avg_n_Wrap(corr(:,{0:20},{250:280}), (/1,2/))

   ; prob_mam = rtest(corr, dimsizes(q_col(:,1,1)), 0)

   ; copy_VarCoords(q_col_annual(1,:,:), prob_mam)
;;;;


;;;;corr with SST
    sst_annual_3d = conform_dims(dimsizes(rh_annual), sst_annual, (/0,2,3/))
    printVarSummary(sst_annual_3d)

    corr = escorc_n(rh_annual, sst_annual_3d,0, 0)
    copy_VarCoords(q_annual(1,:,:,:), corr)

    corr_mean = dim_avg_n_Wrap(corr(:,{0:10},{250:280}), (/1,2/))


 ;==================================
; wks   = gsn_open_wks("pdf", "每层相对湿度与垂直速度的关系")

wks   = gsn_open_wks("pdf", "每层相对湿度与SST的关系")
;;**************************************************************
;画图
;;***************************************************************

    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False

    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True

    resxy@xyLineColors       :=  (/"black"/)          ; change line color
    ; resxy@xylineColors :=  (/"red","blue","green","black"/)
    resxy@xyMarkLineModes    :=  (/"MarkLines"/)

    resxy@xyDashPatterns := (/0/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16/)
    resxy@xyMarkerSizes = (/0.01/)
    ; resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)
    resxy@xyMarkerColors :=  (/"black"/)

    resxy@gsnLeftStringFontHeightF = 0.02
    ; resxy@gsnLeftString = "Corr (NEP omega and Rh)"
    resxy@gsnLeftString = "Corr (NEP SST and Rh)"

    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    resxy@tmXBMode  =  "Explicit"
    ; resxy@trXMaxF =  1.0
    ; resxy@trXMinF = 0.
    ; resxy@tmXBValues = ispan(0,10,2)*0.1
    ; resxy@tmXBLabels = ispan(0,10,2)*0.1
    ; resxy@tmXBMinorOn = True
    ; resxy@tmXBMinorValues = ispan(0,10,1)*0.1

    resxy@trXMaxF =  0.8
    resxy@trXMinF = -0.8
    resxy@tmXBValues = ispan(-8,8,2)*0.1
    resxy@tmXBLabels = ispan(-8,8,2)*0.1
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(-8,8,1)*0.1



    resxy@tiYAxisString = ""
   
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0.0

    plot = gsn_csm_xy(wks,corr_mean ,plev ,resxy)
    
    ;;添加legend 
    gres = True
    gres@YPosPercent = 10.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 70    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 10.
    lineres@lgDashIndexes   = resxy@xyDashPatterns
    lineres@lgMonoMarkerIndex = False 
    lineres@lgMarkerIndexes = resxy@xyMarkers
    lineres@lgItemTypes =  resxy@xyMarkLineModes
    lineres@lgMarkerIndex = resxy@xyMarkers
                        ; expressed as %, 0->100, length of line
    textres = True
    textres@lgLabels = (/"correlation"/)
    textres@lgLabelFontHeightF = 0.02
    plot = simple_legend(wks,plot,gres,lineres,textres)
   

draw(plot)
frame(wks)



end 
