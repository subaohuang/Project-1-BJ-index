;#############################################
;calculate evaporation zonal panel 
;#################################################

begin
fils=systemfunc("ls /mnt/i/oaflux/evapr_oaflux_"+"*"+".nc")
print(fils)
f2 = addfiles(fils, "r")

u_r = f2[:]->evapr(:,{-30:30},:)


u_r@_FillValue = toshort(32766)
u_r&lat@units = "degree_north"
u_r&lon@units = "degree_east"
printVarSummary(u_r)

u =  u_r / 365.0
copy_VarMeta(u_r, u)
delete(u_r)
; u_r = where(condtnl_expr, true_value, false_value)
printVarSummary(u)
;;;处理数据
    ;;;
    u_r_an = month_to_annual(u(:,{-5:5},:), 1)
    printVarSummary(u_r_an)


    u_av = dim_avg_n_Wrap(u_r_an, 1)
    printVarSummary(u_av)

    u_av_pi = lonPivot(u_av, 30)


    u_av_1 = dim_avg_n_Wrap(u_av_pi(0:15,:), 0)

    u_av_2 = dim_avg_n_Wrap(u_av_pi(16:,:), 0)


    printVarSummary(u_av_1)
    ;;;

    u_clm = dim_avg_n_Wrap(u , 0)


;;;

;;;draw plot 
    ;;;; draw1 plane 
    wks = gsn_open_wks("pdf", "/mnt/f/share/era/data/eva_along—equator")
    res = True
    res@xyLineColors=(/"green","blue"/)
    res@xyLineThicknesses=(/3,3/)
    res@xyDashPatterns=(/0,1/)

    res@tmXBMode="Explicit"
    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
              ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7
    plot = gsn_csm_xy(wks, u_av&lon, (/u_av_1,u_av_2/), res)
    ;;;

    ;;;draw climatology plane 
    wks1 = gsn_open_wks("pdf", "/mnt/f/share/era/data/eva—spatial pattern")

    res2 = True
    res2@gsnFrame=False
    res2@gsnDraw=False

    res2@cnFillOn=True
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res2@cnMinLevelValF       =   0
    res2@cnMaxLevelValF       =   5
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

    plot1 = gsn_csm_contour_map(wks1, u_clm, res2)

    draw(plot1)
    frame(wks1)
end