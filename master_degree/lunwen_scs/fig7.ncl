function  one_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end


;;;********************************************************************
;;;主程序
;;;******************************************************************
begin


xi_model=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
dong_model=(/"MIROC5","MPI-ESM-MR"/)
; dong_model=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
; dong_model=(/"MIROC5"/)
;;;全球
latmin=-50
latmax=50
lonmin=0
lonmax=360
lat_sst=fspan(49.5, -49.5, 100)
lon_sst=fspan(0.5, 359.5, 360)
;;;

sst_ximodel=new((/dimsizes(xi_model),30,dimsizes(lat_sst),dimsizes(lon_sst)/),float)
sst_dongmodel=new((/dimsizes(dong_model),30,dimsizes(lat_sst),dimsizes(lon_sst)/),float)


 do i=0,1
    if i.eq.0 then 
     dim=dimsizes(xi_model)   
     models=xi_model
 else
    dim=dimsizes(dong_model)
    models=dong_model
   end if

 do j=0,dim-1
    
    model=models(j)
    yearstart=197901
    yearend=200812
    
    year1=1979
    year2=2008
    
   
    print("******************"+model+"*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_rcp45_r1i1p1_"+"*"+".nc")

    printVarSummary(fils1)
    printVarSummary(fils2)

    printVarSummary(fils1)
    printVarSummary(fils2)
    print(fils1)
    print(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "join")
    delete(fils1)
    delete(fils2)

    vNames= getfilevarnames(f1[0]) ; get variable names of all groups on file
    print(vNames)

    if (any(vNames.eq."LONN179_180")) then 
    time=f1[:]->TIME
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    
  
    sst_r=lonFlip(f1[:]->TOS(timeind,:,:))
    sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r)
    dNames = getvardims(f1[:]->TOS)
      delete(f1)
    print(dNames)
    sst_r!0="time"
    sst_r!1="lat"
    sst_r!2="lon"
    end if

  
    delete(vNames)
    delete(dNames)

    delete(sst_r&time)
    sst_r&time=date(timeind)

    
    sst_mean_annual=month_to_annual(sst_r(:,{latmax:latmin},{lonmin:lonmax}), 1)
    
    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)

    delete(sst_r)
    delete(time)
    delete(date)
    sst=sst_mean_annual
    delete(sst_mean_annual)
    printVarSummary(sst)
    delete(timeind)
    ;;;;

;;;***********************************************************************************************
;;;
;;;;'**********************************************************************************************'

    if i.eq.0 then 
    sst_ximodel(j,:,:,:)=sst

   else
    sst_dongmodel(j,:,:,:)=sst
    end if
    delete([/sst/])
    end do
    delete(dim)
    delete(models)
    end do


    sst_xi_mean=dim_avg_n(sst_ximodel,0)

    sst_dong_mean=dim_avg_n(sst_dongmodel,0)

    printVarSummary(sst_ximodel)
    copy_VarMeta(sst_ximodel(1,:,:,:), sst_xi_mean)
    sst_xi_mean&time=ispan(year1, year2, 1)
    copy_VarCoords(sst_xi_mean, sst_dong_mean)
    printVarSummary(sst_xi_mean)




    print("*****************计算趋势************************")  

    sst_rc1_xi_mean=regCoef_n(ispan(1,dimsizes(sst_xi_mean&time),1),sst_xi_mean(:,:,:),0,0)  
    copy_VarCoords(sst_xi_mean(0,:,:), sst_rc1_xi_mean)
    ; delete(sst_xi_mean)
    printVarSummary(sst_rc1_xi_mean)

    sst_rc1_dong_mean=regCoef_n(ispan(1,dimsizes(sst_dong_mean&time),1),sst_dong_mean(:,:,:),0,0)  
    copy_VarCoords(sst_dong_mean(0,:,:), sst_rc1_dong_mean)
    ; delete(sst_dong_mean)
    printVarSummary(sst_rc1_dong_mean)
 
   

    prob1=one_test(sst_rc1_xi_mean)
    prob2=one_test(sst_rc1_dong_mean)
 
    sst_rc1_xi_mean=(/sst_rc1_xi_mean/)*10.0
    sst_rc1_dong_mean=(/sst_rc1_dong_mean/)*10.0

    
    ;;
    printVarSummary(sst_dongmodel)

 
    sst_rc1_dong1=regCoef_n(ispan(1,dimsizes(sst_dong_mean&time),1),sst_dongmodel(0,:,:,:),0,0)  
    sst_rc1_dong2=regCoef_n(ispan(1,dimsizes(sst_dong_mean&time),1),sst_dongmodel(1,:,:,:),0,0)    
    sst_rc1_dong_mean1=(/(sst_rc1_dong1+sst_rc1_dong2)/2/)
    sst_rc1_dong_mean1=(/sst_rc1_dong_mean1/)*10.0
    copy_VarCoords(sst_rc1_dong_mean, sst_rc1_dong_mean1)
    delete(sst_xi_mean)
    delete(sst_dong_mean)






;;*********************************************************************************************
;;读入观测的海温
;;;*********************************************************************************************
    yearstart=197901
    yearend=200812

;;;
    print("*****************读入海温场*********************")
   ;;;海温场
    fils1=systemfunc("ls /mnt/f/share/ersst/sst.mnmean.v4.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre=f1->sst(timeind,{-50:50},:)

    delete(sst_mre&time)
    printVarSummary(sst_mre)
    sst_mre&time =  date(timeind)
    delete(date)
    delete(f1)
    delete(fils1)
    sst_mre_annual=month_to_annual(sst_mre, 1)
    sst_mre_annual!0="time"
    sst_mre_annual&time=ispan(1979,2008,1)
    printVarSummary(sst_mre_annual)
    

 
  
    ;;;
    print("*****************计算趋势************************")  
    sst_rc1_mre=regCoef_n(ispan(1,dimsizes(sst_mre_annual&time),1),sst_mre_annual(:,:,:),0,0)  
    copy_VarCoords(sst_mre_annual(0,:,:), sst_rc1_mre)
    printVarSummary(sst_rc1_mre)

    prob1_mre=one_test(sst_rc1_mre)

    sst_rc1_mre=(/sst_rc1_mre/)*10

    ;;;;读入hadisst 的海温

     print("*****************读入海温场*********************")
   ;;;海温场
    fils1=systemfunc("ls /mnt/f/share/ersst/HadIsst_sst.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart) 
    sst_mre_had=lonFlip(f1->sst(timeind,{-50:50},:))
    delete(sst_mre_had&time)
    printVarSummary(sst_mre)
    sst_mre_had&time = date(timeind)
    delete(date)
    delete(f1)
    delete(fils1)
    sst_mre_had_annual=month_to_annual(sst_mre_had, 1)
    sst_mre_had_annual!0="time"
    sst_mre_had_annual&time=ispan(1979,2008,1)
    printVarSummary(sst_mre_had_annual)

    sst_cha= dim_avg_n_Wrap(sst_mre_had_annual(:,{-5:5},{80:160}),(/1,2/)) - dim_avg_n_Wrap(sst_mre_had_annual(:,{-5:5},{200:280}),(/1,2/))
 


    rc_sst=regCoef_n(ispan(1,30,1),sst_cha , 0, 0)
    rc_sst=(/rc_sst*10/)
    print("rc_sst="+rc_sst)
   
 
  
    ;;;
    print("*****************计算趋势************************")  
    sst_rc1_had=regCoef_n(ispan(1,dimsizes(sst_mre_had_annual&time),1),sst_mre_had_annual(:,:,:),0,0)  
    copy_VarCoords(sst_mre_had_annual(0,:,:), sst_rc1_had)
    printVarSummary(sst_rc1_had)

    prob1_had=one_test(sst_rc1_had)

    sst_rc1_had=(/sst_rc1_had/)*10











    print("**********************画图*****************************")
    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/scs_fig.7")


    colormap=read_colormap_file("GMT_jet")
    printVarSummary(colormap)
    cmap=new((/14,4/),typeof(colormap))
    printVarSummary(colormap)

    cmap(0:6,:)=colormap(2:92:15,:)
    ; cmap(0:6,:)=colormap(2:60:8,:)
    cmap(7:13,:)=colormap(130:250:20,:)
    

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False



    res@gsnCenterStringFontHeightF=0.02
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False


    ; res@cnFillPalette="BlueRed"
    res@cnFillPalette=cmap

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnFillDrawOrder="PreDraw"
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05

    
    
    res@tmXTOn=True
    res@tmYROn=True
    res@tmXUseBottom=False
    res@tmXBLabelsOn=False
    res@tmXTMinorOn=False
    res@tmXBMinorOn=False
    res@tmYLMinorOn=False

    

    res@tmXBLabelFontHeightF=0.012
    res@tmYLLabelFontHeightF=0.012

    res@tmXBMajorOutwardLengthF=0
    res@tmXBMajorLengthF=0.007

    res@tmXTMajorOutwardLengthF=0
    res@tmXTMajorLengthF=0.007

    res@tmYLMajorOutwardLengthF=0
    res@tmYLMajorLengthF=0.01

    res@tmYLMajorThicknessF=2.0
    res@tmXTMajorThicknessF=2.0
    res@tmBorderThicknessF=2.5

    
    res@tmXTMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXTValues=ispan(0,360,60)
    res@tmXTLabelsOn=True
    res@tmXTLabels=(/"0E","60~S~o~N~E","120~S~o~N~E","180~S~o~N~E","120~S~o~N~W","60~S~o~N~W","0E"/)
    res@tmXTLabelDeltaF=-0.6
    
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXBValues=ispan(0,360,60)
    

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLValues=ispan(50,-50,25)
    res@tmYLLabels=(/"50~S~o~N~N","25~S~o","0","20~S~o","50~S~o~N~S"/)
    res@tmYLLabelDeltaF=-0.6


  
    res@mpFillDrawOrder="PostDraw"
    res@mpLandFillColor="gray"
    res@mpFillOn=True
    
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = latmin
    res@mpMaxLatF             = latmax
    res@mpMinLonF             = lonmin
    res@mpMaxLonF             = lonmax
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpShapeMode="FreeAspect"
    res@vpWidthF=0.65
    res@vpHeightF=0.35

    plot=new(4,graphic)
    plot(0) = gsn_csm_contour_map(wks1,sst_rc1_mre,res)
    res@tmYLLabelsOn=False
    plot(1) = gsn_csm_contour_map(wks1,sst_rc1_had,res)
    res@tmYLLabelsOn=True
    res@tmXTLabelsOn=False
    plot(2) = gsn_csm_contour_map(wks1,sst_rc1_xi_mean,res)
    res@tmXTLabelsOn=False
    res@tmYLLabelsOn=False
    plot(3) = gsn_csm_contour_map(wks1,sst_rc1_dong_mean1,res)
     
    
         ;;;打点
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False

        iplot=new(4,graphic)

        iplot(2)=gsn_csm_contour(wks1,prob1,res_pro)
        iplot(3)=gsn_csm_contour(wks1,prob2,res_pro)
        iplot(1)=gsn_csm_contour(wks1,prob1_had,res_pro)
        iplot(0)=gsn_csm_contour(wks1,prob1_mre,res_pro)
        
        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 1.0         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0025      ; make dots larger
        do i=0,3
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
        ;;;  
        overlay(plot(i), iplot(i))
        end do

       ;;;添加N与S，热带海洋学报要求
    
    ;;位置0.784909 ，0.7849-0.247587 ， 0.5096, 0.50-0.24758
    ;字体高度0.006


    ; txres=True
    ; txres@txFontHeightF=0.008
    ; gsn_text_ndc(wks1, "N", 0.011, 0.784909-0.01, txres)
    ; gsn_text_ndc(wks1, "S", 0.011, 0.784909-0.247587+0.01, txres)
    
    ; gsn_text_ndc(wks1, "N", 0.011, 0.5096-0.01, txres)
    ; gsn_text_ndc(wks1, "S", 0.011, 0.5096-0.247587+0.01, txres)


; drawNDCGrid(wks1)
        pres=True
        pres@gsnPanelLabelBar=True
        pres@pmLabelBarWidthF = 0.1
        pres@pmLabelBarHeightF = 0.45
        pres@lbLabelStride = 2
        pres@lbOrientation = "Vertical"
        pres@gsnPanelDebug=True

        pres@gsnPanelFigureStrings=(/"a","b","c","d"/)
        pres@gsnPanelFigureStringsJust="TopLeft"

        txres=True
        txres@txFontHeightF = 0.01
        lbTitleString = "(k~F34~W~F0~decade~S~-1~N~)"
        gsn_text_ndc(wks1,lbTitleString,.95,0.73,txres)

        gsn_panel(wks1, plot, (/2,2/), pres)


end
    
  
