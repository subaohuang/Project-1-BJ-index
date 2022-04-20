;#############################################
;calculate u850 zonal panel 
;#################################################

f2=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
date=cd_calendar(f2->time, 1)
timeind=ind(date.ge.197901.and.date.le.201312)
u_r = f2->precip(timeind,{-5:5},:)


u_r_an = month_to_annual( u_r, 1)
printVarSummary(u_r_an)

u_av = dim_avg_n_Wrap(u_r_an, 1)
printVarSummary(u_av)

u_av_pi = lonPivot(u_av, 30)


u_av_1 = dim_avg_n_Wrap(u_av_pi(0:15,:), 0)

u_av_2 = dim_avg_n_Wrap(u_av_pi(16:,:), 0)


printVarSummary(u_av_1)

;;;; draw 
 wks = gsn_open_wks("pdf", "/mnt/f/share/era/data/pre_along—equator")
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
;;

plot = gsn_csm_xy(wks, u_av&lon, (/u_av_1,u_av_2/), res)
