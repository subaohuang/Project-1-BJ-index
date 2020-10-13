;#################################################
;calcprelate evaporation  - precipitation zonal panel 
;#################################################

begin

;;read  eva data 
    fils=systemfunc("ls /mnt/i/oaflux/evapr_oaflux_"+"*"+".nc")

    f2 = addfiles(fils, "r")
    u_r = f2[:]->evapr(:,:,:)

    u_r@_FillValue = toshort(32766)
    u_r&lat@units = "degree_north"
    u_r&lon@units = "degree_east"
    printVarSummary(u_r)

    u =  u_r / 365.0
    copy_VarMeta(u_r, u)
    delete(u_r)

    printVarSummary(u)
    



    ;;;处理数据
        ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

        yo = fspan(-88.75, 88.75, 72)
        xo = fspan(1.25, 358.75, 144)

        u_new = linint2_Wrap(u&lon, u&lat, u, True, xo, yo, 0)

        delete(u_new&LAT)
        delete(u_new&LON)
        u_new!1 = "lat"
        u_new!2 = "lon"
        u_new&lat = yo
        u_new&lon = xo
        printVarSummary(u_new)
     

        ;;;

       
        u_av = dim_avg_n_Wrap(u_new(:,{-5:5},:), 1)
  
    
        u_av_pi = lonPivot(u_av, 30)
    
       
        u_av_1 = dim_avg_n_Wrap(u_av_pi(0:15*12-1,:), 0)

        u_av_2 = dim_avg_n_Wrap(u_av_pi(15*12:,:), 0)
       
        ;;116.25存在插值后缺侧的现象，现在插值下
        m =  linmsg(u_av_1({100:120}),0)
        copy_VarCoords(u_av_1({100:120}), m)
        u_av_1({116.25}) = m({116.25})
        delete(m)

        m =  linmsg(u_av_2({100:120}),0)
        copy_VarCoords(u_av_2({100:120}), m)
        u_av_2({116.25}) = m({116.25})



        ;;;
        u_clm = dim_avg_n_Wrap(u_new(:,{-30:30},:) , 0)
    ;;;
;;


;;;read pre data

    f1=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    pre_r = f1->precip(timeind,:,:)


    pre_r_an = month_to_annual( pre_r(:,{-5:5},:), 1)

    pre_av = dim_avg_n_Wrap(pre_r_an, 1)

    pre_av_pi = lonPivot(pre_av, 30)

    pre_av_1 = dim_avg_n_Wrap(pre_av_pi(0:15,:), 0)

    pre_av_2 = dim_avg_n_Wrap(pre_av_pi(16:,:), 0)

    pre_clm = dim_avg_n_Wrap(pre_r(:,{-30:30},:), 0)

;;;

;;;calculate E-P

   ep1 = pre_av_1 - u_av_1 
   ep2 = pre_av_2 - u_av_2 
   copy_VarCoords(pre_av_1, ep1)
   copy_VarCoords(pre_av_2, ep2)

   ep_clm = pre_clm - u_clm
   copy_VarCoords(pre_clm, ep_clm)
;;;


;;;

;;;draw plot 
    ;;;; draw1 plane 
    wks = gsn_open_wks("pdf", "/mnt/f/share/era/data/ep_along—eqpreator")
    res = True
    res@xyLineColors=(/"green","blue"/)
    res@xyLineThicknesses=(/3,3/)
    res@xyDashPatterns=(/0,1/)
     
    res@trXMinF = 28.75
    res@trXMaxF = 386.25
    res@tmXBMode="Explicit"
    res@tmXBValues=(/30,90,150,210,270,330,386.25/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
              ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7

    plot = gsn_csm_xy(wks, u_av_1&lon, (/ep1,ep2/), res)
    ;;;

    ;;;draw climatology plane 
    wks1 = gsn_open_wks("pdf", "/mnt/f/share/era/data/ep—spatial pattern")

    res2 = True
    res2@gsnFrame=False
    res2@gsnDraw=False

    res2@cnFillOn=True
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode = "ManualLevels"      ; set manpreal contoprer level
    res2@cnMinLevelValF       =   0
    res2@cnMaxLevelValF       =   8
    res2@cnLevelSpacingF      =  0.5
    res2@cnLineLabelsOn=False
    res2@lbLabelBarOn=True

    res2@mpLimitMode           = "LatLon"
    res2@mpMinLatF             = -30
    res2@mpMaxLatF             = 30
    res2@mpMinLonF             =0
    res2@mpMaxLonF             = 360
    res2@mpCenterLonF          = (res2@mpMinLonF+res2@mpMaxLonF)/2.0

    res2@mpShapeMode="FreeAspect"
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7 

    plot1 = gsn_csm_contour_map(wks1, ep_clm, res2)

    draw(plot1)
    frame(wks1)
end