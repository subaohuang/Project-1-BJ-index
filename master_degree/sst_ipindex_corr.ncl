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
    diri="/mnt/f/share/era/data/"
    ipindex=asciiread("/mnt/f/share/era/data/ip_index_1950_2010.ascii", 732, "float")

    ipindex_annual=month_to_annual(ipindex, 1)
    ipindex_anu_std=dim_standardize_n(ipindex_annual,0,0)
    ipindex_anu_std!0="time"
    ipindex_anu_std&time=ispan(1,61,1)
     
    ipindex_anu_std_de=dtrend(ipindex_anu_std, True)
    copy_VarCoords(ipindex_anu_std, ipindex_anu_std_de)
    
    ipindex_anu_std_sm=runave_n_Wrap(ipindex_anu_std, 9, 0, 0)

    fils2=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc") 
    f2=addfile(fils2, "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201012.and.date.ge.195001)
    sst=lonFlip(f2->sst(timeind,{60:-40},:))
    delete(sst&time)

    sst&time=date(timeind)
    
    sst_annual=month_to_annual(sst, 1)
    sst_annual&year=ispan(1, 61, 1)

    sst_de=dtrend_msg_n(ispan(0,dimsizes(sst_annual&year)-1,1), sst_annual, True, False, 0)
    copy_VarCoords(sst_annual, sst_de)
  
   
    sst_sm=runave_n_Wrap(sst_annual, 9, 0, 0)
    
    r=escorc_n(sst_sm, ipindex_anu_std_sm, 0, 0)

   
    copy_VarCoords(sst_sm(0,:,:), r)
    
    printVarSummary(r)

     ; n1=equiv_sample_size(sst_sm(latitude|:,longitude|:,year|:), 0.05, 0)
     ; n2=equiv_sample_size(ipindex_anu_std_sm, 0.05, 0)

  ;;;做显著性检验 计算有效自由度
    r1=esacr(sst_sm(latitude|:,longitude|:,year|:), 2)  ;计算自相关系数
    r2=esacr(ipindex_anu_std_sm, 2)  ;计算自相关系数
   
    n_o=61 
    n1=(/n_o*(1-r1(:,:,1)*r1(:,:,2))/(1+r1(:,:,1)*r1(:,:,2))/)
    n2=61*(1-r2(1)*r2(2))/(1+r2(1)*r2(2))  
    
    n2nd=conform_dims(dimsizes(n1), n2, 0)
    n=(n1+n2nd)/2
    n_in=floattointeger(n)
   
    prob= rtest(r, n_in+2, 0)
    copy_VarCoords(r, prob)
  ;;;

    
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/sst&指数的9年平滑的相关性")
    gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
    res@gsnLeftString="w trend"
    res@gsnRightString="1950-2010"

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnFillMode="RasterFill"
    res@cnRasterSmoothingOn=True
    

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.8
    res@cnMaxLevelValF       =   0.8
    res@cnLevelSpacingF      =  0.1

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(r&latitude)
    res@mpMaxLatF             = max(r&latitude)
    res@mpMinLonF             = min(r&longitude)
    res@mpMaxLonF             = max(r&longitude)
    res@mpCenterLonF          =avg(r&longitude)

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.8
    res@vpHeightF = 0.2

    ; res@cnMinLevelValF       =   -5.0
    ; res@cnMaxLevelValF       =   5.0
    ; res@cnLevelSpacingF      =  1.0

    plot= gsn_csm_contour_map(wks,r,res)



    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05      ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot=gsn_csm_contour(wks,prob,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger


    iplot = gsn_contour_shade(iplot, 0,0.05, opt)

    ;;;  

    overlay(plot, iplot)

    pres=True
    pres@gsnPanelLabelBar=True

    draw(plot)
    frame(wks)
end