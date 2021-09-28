function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1)/),float)
    do i=0,11   
    data_12(i,:,:) = data(i::12,:)
    end do 
return data_12
end 

begin
    yearstart = 25601
    yearend = 28012

    lat1 = 0
    lat2 = 20
    lon1 = 250
    lon2 = 280

;;;;read contrl data
    modelname = "B2000_alt_north_year_CTRL"
    diri = "/home/yangsong3/data-observation/linshh/CESM-data/B/"+modelname+"/"
    filename= modelname+".cam.h1.0251-0280_chazhi.nc"

    f1 = addfile(diri+filename, "r")
   
    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)   

    plev = f1->lev({1000:300})


    t_ct = f1->T(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q_ct = f1->Q(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    ; q = q/1000.0  ;;;convert to kg/kg
    hgt_ct = f1->Z3(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;


;;;;read data
    diri = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year/"
    filename="B2000_alt_north_year.cam.h1.0251-0280_chazhi.nc"

    f1 = addfile(diri+filename, "r")
   
    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)

    t_sen = f1->T(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q_sen = f1->Q(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    ; q = q/1000.0  ;;;convert to kg/kg
    hgt_sen = f1->Z3(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;

    Cp = 1004
    Lv = 2.5*10^6
    g = 9.8
    
    hgt_ct = (/hgt_ct*9.8/10^3/)
    t_ct  = Cp*t_ct/10^3
    q_ct  = Lv*q_ct/10^3

    hgt_sen = (/hgt_sen*9.8/10^3/)
    t_sen  = Cp*t_sen/10^3
    q_sen  = Lv*q_sen/10^3


    season  = (/"MAM"/)
  

;;;;;season 

    t_ct_nep = dim_avg_n_Wrap(t_ct, (/2,3/))
    q_ct_nep = dim_avg_n_Wrap(q_ct, (/2,3/))
    hgt_ct_nep = dim_avg_n_Wrap(hgt_ct, (/2,3/))

    t_ct_12 = split_12(t_ct_nep)
    q_ct_12 = split_12(q_ct_nep)
    hgt_ct_12 = split_12(hgt_ct_nep)

    t_ct_nep_MAM  = dim_avg_n_Wrap(t_ct_12(0:4,:,:),(/0/))
    q_ct_nep_MAM  = dim_avg_n_Wrap(q_ct_12(0:4,:,:),(/0/))
    hgt_ct_nep_MAM  = dim_avg_n_Wrap(hgt_ct_12(0:4,:,:),(/0/))

    t_sen_nep = dim_avg_n_Wrap(t_sen, (/2,3/))
    q_sen_nep = dim_avg_n_Wrap(q_sen, (/2,3/))
    hgt_sen_nep = dim_avg_n_Wrap(hgt_sen, (/2,3/))

    t_sen_12 = split_12(t_sen_nep)
    q_sen_12 = split_12(q_sen_nep)
    hgt_sen_12 = split_12(hgt_sen_nep)

    t_sen_nep_MAM  = dim_avg_n_Wrap(t_sen_12(0:4,:,:),(/0/))
    q_sen_nep_MAM  = dim_avg_n_Wrap(q_sen_12(0:4,:,:),(/0/))
    hgt_sen_nep_MAM  = dim_avg_n_Wrap(hgt_sen_12(0:4,:,:),(/0/))
    
    t_ct_nep_clm_mam = dim_avg_n_Wrap(t_ct_nep_MAM, 0)
    q_ct_nep_clm_mam = dim_avg_n_Wrap(q_ct_nep_MAM, 0)
    hgt_ct_nep_clm_mam = dim_avg_n_Wrap(hgt_ct_nep_MAM, 0)
    
    t_sen_nep_clm_mam = dim_avg_n_Wrap(t_sen_nep_MAM, 0)
    q_sen_nep_clm_mam = dim_avg_n_Wrap(q_sen_nep_MAM, 0)
    hgt_sen_nep_clm_mam = dim_avg_n_Wrap(hgt_sen_nep_MAM, 0)
    

    t_nep_change = t_ct_nep_clm_mam
    t_nep_change = t_sen_nep_clm_mam - t_ct_nep_clm_mam

    q_nep_change = q_ct_nep_clm_mam
    q_nep_change = q_sen_nep_clm_mam - q_ct_nep_clm_mam

    hgt_nep_change = hgt_ct_nep_clm_mam
    hgt_nep_change = hgt_sen_nep_clm_mam - hgt_ct_nep_clm_mam
   

    mse_change = hgt_nep_change + t_nep_change + q_nep_change

    print("mse_change="+sum(mse_change))
    
    ; printMinMax(mse_nep_clm, 1)
    ; exit()

    ; order = (/"(c)","(d)"/)

    plot = new(1,graphic)
    wks = gsn_open_wks("pdf", "B2000-ATL-YEAR-MSE-CPT-Lvq_change")
;;;;
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       :=  (/"red","blue","green","black"/)          ; change line color
    ; resxy@xylineColors :=  (/"red","blue","green","black"/)
    resxy@xyMarkLineModes    :=  (/"MarkLines","MarkLines","MarkLines","MarkLines"/)
    resxy@xyDashPatterns := (/0,0,0,0/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16,16,16,16/)
    resxy@xyMarkerSizes = (/0.01,0.01,0.01,0.01/)
    ; resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)
    resxy@xyMarkerColors :=  (/"red","blue","green","black"/)


    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    resxy@tmXBMode  =  "Explicit"
 
    resxy@tiYAxisString = ""
   
    resxy@trYReverse        = True                       ; reverse Y-axis

    resxy@trXMaxF :=  2.
    resxy@trXMinF := -1.
    resxy@tmXBValues := ispan(-10,20,5)*0.1
    resxy@tmXBLabels := ispan(-10,20,5)*0.1
    resxy@tmXBMinorOn = False
    ; resxy@tmXBMinorValues := ispan(-10,40,1)*0.1
    resxy@tiXAxisString = "(10~S~3~N~J/kg)"
    resxy@gsnLeftString = "(b) MSE change profile"

    plot = gsn_csm_xy(wks, (/t_nep_change,q_nep_change, \
                             hgt_nep_change,mse_change/),plev ,resxy)
    
 
    ;;添加legend 
    gres = True
    gres@YPosPercent = 20.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 80    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = resxy@xyDashPatterns
    lineres@lgMonoMarkerIndex = False 
    lineres@lgMarkerIndexes = resxy@xyMarkers
    lineres@lgItemTypes =  resxy@xyMarkLineModes
    lineres@lgMarkerIndex = resxy@xyMarkers
                        ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"CpT","Lq","hgt","MSE"/)
    textres@lgLabelFontHeightF = 0.01
    plot = simple_legend(wks,plot,gres,lineres,textres) 
    


    draw(plot)
    frame(wks)

    ; pres = True
    ; gsn_panel(wks, plot, (/2,2/), pres)

end




