   
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
   diri="/mnt/f/share/NOAA/"
   fils1=systemfunc("ls "+diri+"olr.mon.mean.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    year1=1979
    year2=2012
    yearend=201212
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    olr_r=f1->olr(timeind,:,:)
    olr=(/olr_r@add_offset+olr_r@scale_factor*olr_r/)
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    delete(olr&time)
    olr&time=date(timeind)
    printVarSummary(olr)
    olr_pi=lonPivot(olr, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "DJF")
    olr_djf=olr_djf_r(0:,:,:)
    
    olr_djf_pool=dim_avg_n_Wrap(olr_djf(:,{10:10},{80:150}), (/1,2/))
    olr_djf_pool=(/dim_standardize_n(olr_djf_pool, 0, 0)*-1.0/)
    
    ;;;
    delete(date)
    delete(timeind)

    f2=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/mnt/f/share/GPCC/precip.mon.total.v7.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi=lonPivot(p, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    p_djf_r=month_to_season(p_pi(time|:,lat|:,lon|:), "DJF")
    ;;注意点
    p_djf=p_djf_r(0:,:,:)
    
    p_djf_pool=dim_avg_n_Wrap(p_djf(:,{10:10},{80:150}), (/1,2/))
    p_djf_pool=(/dim_standardize_n(p_djf_pool, 0, 0)/)
    delete(date)
    delete(timeind)
    delete(p)
    delete(p_pi)
    ;;;
    
    ;;;GPCC 陆地降水
    f3=addfile("/mnt/f/share/GPCC/precip.mon.total.v7.nc", "r")
    date=cd_calendar(f3->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    p=f3->precip(timeind,:,:)
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi=lonPivot(p, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    p_djf_r_land=month_to_season(p_pi(time|:,lat|:,lon|:), "DJF")
    ;;注意点
    p_djf_land=p_djf_r_land(0:,:,:)
    printVarSummary(p_djf_land)
    
    p_djf_land_pool=dim_avg_n_Wrap(p_djf_land(:,{10:10},{80:150}), (/1,2/))
    p_djf_land_pool=(/dim_standardize_n(p_djf_land_pool, 0, 0)/)
     delete(date)
    delete(timeind)
    ;;;

    ;;;Icoads 云量
    f4=addfile("/mnt/f/share/ICOADS/cldc.mean.nc", "r")
    date=cd_calendar(f4->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    cldc_r=f4->cldc(timeind,:,:)
    cldc=cldc_r*cldc_r@scale_factor+cldc_r@add_offset
    copy_VarCoords(cldc_r, cldc)
    delete(cldc&time)
    cldc&time=date(timeind)
    printVarSummary(cldc)
    
    ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    opt=(/0,5/)
    ;;先插值lon方向的
    cldc=(/linmsg_n(cldc, opt, 2)/)
    ;;插值lat方向的
    opt=(/0,1/)
    cldc=(/linmsg_n(cldc, opt, 1)/)
    delete(opt)

    cldc_pi=lonPivot(cldc, 30)
    cldc_djf_r=month_to_season(cldc_pi(time|:,lat|:,lon|:), "DJF")
    ;;注意点
    cldc_djf=cldc_djf_r(0:,:,:)
   
    cldc_djf_mean=dim_avg_n_Wrap(cldc_djf(:,{30:-30},:), (/1,2/))
    do i=0,dimsizes(cldc_djf&time)-1
        cldc_djf(i,:,:)=(/cldc_djf(i,:,:)-cldc_djf_mean(i)/)
    end do


    cldc_djf_pool=dim_avg_n_Wrap(cldc_djf(:,{10:10},{80:150}), (/1,2/))
    cldc_djf_pool=(/dim_standardize_n(cldc_djf_pool, 0, 0)/)
    ;;;;
   
    ;;;three year runmean

    olr_djf_pool=(/runave(olr_djf_pool, 3, 0)/)
    p_djf_pool=(/runave(p_djf_pool, 3, 0)/)
    p_djf_land_pool=(/runave(p_djf_land_pool, 3, 0)/)
    cldc_djf_pool=(/runave(cldc_djf_pool, 3, 0)/)
    
    p_djf_pool@_FillValue=olr_djf_pool@_FillValue
    p_djf_land_pool@_FillValue=olr_djf_pool@_FillValue


    ;;;
    wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/OLR与GPCP,cloud的时间序列")
;;**************************************************************
;画图
;;***************************************************************


    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnYRefLine=0.0
    
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@trXMaxF=34
    res@trXMinF=1
    res@trYMaxF=3.0
    res@trYMinF=-3.0
    res@tmXBValues=ispan(1,34,10)
    res@tmXBLabels=ispan(1979,2012,10)
    res@tmXBMinorValues=ispan(1,34,1)
    
    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines

    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@xyLineColors=(/"blue","black","black","red"/)
    res@xyDashPatterns=(/0,0,2,0/)
    res@xyLineThicknesses=(/2,2,2,2/)
    res@vpHeightF=0.4
    res@vpWidthF=0.6





    
    plot=gsn_csm_xy(wks, ispan(1, 34, 1), (/olr_djf_pool,p_djf_pool,p_djf_land_pool,cldc_djf_pool/), res)
    

 ;;添加legend 
    gres = True
    gres@YPosPercent = 96.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 6    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    lineres@lgLineColors = (/"blue","black","black","red"/)
    lineres@lgDashIndexes=(/0,0,2,0/)
    lineres@lgLineThicknesses = 2.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"-OLR","GPCP","GPCC","Cloud"/)
    textres@lgLabelColors=(/"blue","black","black","red"/)

    plot = simple_legend(wks,plot,gres,lineres,textres)   

  




    draw(plot)
    frame(wks)
    end