   begin
   diri="/mnt/f/share/ncep/ncep1/"
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
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf=month_to_season(olr(time|:,lat|:,lon|:), "DJF")

    olr_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    olr_djf(0:dimsizes(olr_djf&time)-2,:,:)=(/olr_djf(1:,:,:)/)
    olr_djf(dimsizes(olr_djf&time)-1,:,:)=(/olr({yearend},:,:)/)
    ;;;

    olr_djf1=dim_avg_n_Wrap(olr_djf({1979:1994},{-30:30},:), 0)
    olr_djf2=dim_avg_n_Wrap(olr_djf({1995:},{-30:30},:), 0)
    olr_diff=olr_djf2-olr_djf1

  wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/ncep_olr时间段比较")
    res=True
    

    res@cnFillOn=True
    res@cnLinesOn=False
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             =0
    res@mpMaxLonF             = 360
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    plot=gsn_csm_contour_map(wks, olr_diff, res)
    end 