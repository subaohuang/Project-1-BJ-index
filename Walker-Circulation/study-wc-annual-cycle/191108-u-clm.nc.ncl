;;;************************************************************
;; draw the annual cycle of u10 
;；探究不同经度的10-m
;; linshuheng 2019.10.16
;;;*********************************************************
function annual_cycle(u10_inten)
local u10_inten_month, nmos,nmo
begin
    u10_inten_month = new(12,typeof(u10_inten))
    
    nmos = 12
    ntim = 456
    do nmo=0,nmos-1
    u10_inten_month(nmo) = dim_avg( u10_inten(nmo:ntim-1:nmos) )
    end do
return u10_inten_month
end
;;---------------------------------------
;;read  wind and slp data
;----------------------------------------
begin
    fils = systemfunc("ls /home/ys17-19/lsh/data/era/era_pre_global.nc")
    f3 = addfile(fils, "r")

    time=f3->time
    date=cd_calendar(time, 1)
    tind=ind(date.le.201612.and.date.ge.197901)
    u10 = short2flt(f3->u(tind,{200},:,:))
    
    u10_clm  = dim_avg_n_Wrap(u10, 0)
     
  
;;;===================画图=================================
   wks=gsn_open_wks("pdf", "u-pattern")
;;;
    res            = True                         ; plot mods desired
   
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlueRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -20
    res@cnMaxLevelValF       =   20
    res@cnLevelSpacingF      =  5

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -60.0
    res@mpMaxLatF             = 60.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    plot_sst = gsn_csm_contour_map(wks, u10_clm, res)
   
   draw(plot_sst)
   frame(wks)

 end

    

