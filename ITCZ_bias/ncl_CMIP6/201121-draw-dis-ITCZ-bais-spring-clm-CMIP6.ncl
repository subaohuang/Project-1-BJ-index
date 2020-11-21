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

model_name = (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
   , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
    "GISS-E2-1-G",  "GISS-E2-1-G-CC","HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)
iplot= plot
plot1 = plot

wks = gsn_open_wks("pdf", "CMIP6模式降水偏差分布")


;;;读入GPCP的数据,计算其春季气候态
    f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.200412)
    p=f1->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
  
    p_annual = month_to_season(p, "MAM")
    p_annual = (/dtrend_msg_n(ispan(1,dimsizes(p_annual&time),1),p_annual,False,True,0)/)
    p_clm = dim_avg_n_Wrap(p_annual, 0) 
    printVarSummary(p_clm)
    printMinMax(p_clm, 1)
;;;;


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 197901
    yearend = 200412
    fils1 =  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/pr_chazhi/pr_Amon_"+model+"_historical_r1i1p1_"+"*")

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    pr = lonFlip(f1->pr(timeind,:,:))
    
    pr = (/dtrend_msg_n(ispan(1,dimsizes(pr&time),1),pr,False,True,0)/)
    pr = (/pr*86400/) ;;conver to mm/day
    printVarSummary(pr)
    printMinMax(pr, 1)



    pr_annual = month_to_season(pr, "MAM")

    PR_clm = dim_avg_n_Wrap(pr_annual, 0)

    diff = PR_clm
    diff = (/PR_clm - p_clm/) 
 
   printMinMax(diff, 1)

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
    res@cnMinLevelValF       =   -3.0
    res@cnMaxLevelValF       =   3.0
    res@cnLevelSpacingF      =  0.5
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
    plot(k) = gsn_csm_contour_map(wks, diff,res)

 print("********************************finish calculate"+" "+model+"************************************")
 end do 


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarHeightF = 0.05
pres@pmLabelBarWidthF  = 0.6
pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbOrientation = "Vertical"
pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/6,5/), pres)
end 



