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


begin
;读入数据

    read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
    print(read_name)
    model_name = read_name 

    dim_model = dimsizes(model_name)

    plot = new(dim_model+1,graphic)
    iplot= plot
    plot1 = plot

    yearstart = 197901
    yearend = 200412

    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)

;;;


;;;读入GPCP的数据,计算其春季气候态
    f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.200412)
    p=f1->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)

    p = (/rmMonAnnCycTLL(p)/)
 
    ; p_annual = month_to_season(p, "MAM")

    p_annual = month_to_season(p, "DJF")
    p_annual(0:dimsizes(p_annual&time)-2,:,:) = p_annual(1:,:,:)
    p_annual(dimsizes(p_annual&time)-1,:,:) = p(dimsizes(p&time)-1,:,:)

    
    p_annual = (/dtrend_msg_n(ispan(1,dimsizes(p_annual&time),1),p_annual,False,True,0)/)

    p_clm = dim_avg_n_Wrap(p_annual, 0) 
;;;


; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)

      timeind:=ind(date.le.200512.and.date.ge.197901)

      sst =  lonFlip(f1->sst(timeind,:,:))

      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)

      
      ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

      ; ssta_index = dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
      ssta_index = dim_avg_n_Wrap(ssta(:,{5:-5},{210:270}), (/1,2/))
      printVarSummary(ssta_index)

      obs_nino34_r = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      
      obs_nino34 = obs_nino34_r(1:)
      printVarSummary(obs_nino34)
      
      obs_nino34_std = dim_stddev_n_Wrap(obs_nino34, 0)

      ely=(/1980,1983,1987,1988,1992,1995,1998,2003/) - 1
      
      ; obs_nino34_std = dim_stddev_n_Wrap(obs_nino34, 0)
      ; ind_obs_nino = ind(obs_nino34.ge.obs_nino34_std)

      ind_obs_nino = ely - 1979

      obs_nino_intensity = avg(obs_nino34(ind_obs_nino))
    


;;;;;

;;;select enino year spring preciptaion and nino intensity
    printVarSummary(p_annual)
    printVarSummary(obs_nino34)

    obs_p_ninoyear  = dim_avg_n_Wrap(p_annual(ind_obs_nino,:,:) , 0)

;;;
      
;;;


wks = gsn_open_wks("pdf", "30个CMIP5模式中Enino年冬季降水异常分布")


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
   
    fils1 =  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/pr_chazhi/pr_Amon_"+model+"_historical_r1i1p1_"+"*")

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    pr = lonFlip(f1->pr(timeind,:,:))
    pr = (/rmMonAnnCycTLL(pr)/)


    pr = (/dtrend_msg_n(ispan(1,dimsizes(pr&time),1),pr,False,True,0)/)
    
    pr = (/pr*86400/) ;;conver to mm/day
    printVarSummary(pr)
    printMinMax(pr, 1)



    ; pr_annual = month_to_season(pr, "MAM")  ;;;1950年春季

    pr_annual = month_to_season(pr, "DJF")

    pr_annual(0:dimsizes(pr_annual&time)-2,:,:) = pr_annual(1:,:,:)

    pr_annual(dimsizes(pr_annual&time)-1,:,:) = pr(dimsizes(pr&time)-1,:,:)
    PR_clm = dim_avg_n_Wrap(pr_annual, 0)

    ; ; ;;;0 -20 100-150
    ; ; siindex = avg(diff({0:-20},{100:150}))
    ; ; ;;

    ; ;;;0 -20 150-
    ; siindex = avg(diff({0:10},{210:280}))
    ; ;;
    ; dbitcz(k) = siindex 

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    tos = lonFlip(f2[:]->tos(timeind,:,:))
       
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos = (/rmMonAnnCycTLL(tos)/)
    tos_djf = month_to_season(tos, "DJF")

    tos_djf(0:dimsizes(tos&time)/12 -2,:,:) = tos_djf(1:,:,:)
    tos_djf(dimsizes(tos&time)/12 -1,:,:) = tos(dimsizes(tos&time)-1,:,:)
    


    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{190:240}), (/1,2/))
    nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{210:270}), (/1,2/))

    nino3_4_std = dim_stddev_n_Wrap(nino3_4, 0)

    yearind := ind(nino3_4.ge.0.8*nino3_4_std)

    tos_nino := nino3_4(yearind)

    nino_intensity = avg(tos_nino)
    ;;;;
    ; enino_intensity(k) = nino_intensity

    model_p_ninoyear = dim_avg_n_Wrap(pr_annual(yearind,:,:) ,0)


    diff_p = model_p_ninoyear
    
    ; diff_p = (/model_p_ninoyear - obs_p_ninoyear/)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn = False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -4
    res@cnMaxLevelValF       =   4
    res@cnLevelSpacingF      =  1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 120.0
    res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnRightString = ""
    res@gsnLeftString = model
    res@gsnLeftStringFontHeightF = 0.03
    plot(k) = gsn_csm_contour_map(wks, diff_p,res)

 print("********************************finish calculate"+" "+model+"************************************")
end do
    

    plot(dim_model) = gsn_csm_contour_map(wks, obs_p_ninoyear,res)

    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@pmLabelBarHeightF = 0.05
    pres@pmLabelBarWidthF  = 0.6
    pres@pmLabelBarOrthogonalPosF = -0.02
    ; pres@lbOrientation = "Vertical"
    pres@lbLabelFontHeightF = 0.01
    gsn_panel(wks, plot, (/6,6/), pres)
end 


