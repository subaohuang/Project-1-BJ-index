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
   diri="/mnt/f/share/ICOADS/"
    fils1=systemfunc("ls "+diri+"std/cldc.mean.nc")
    
    ; diri="/mnt/f/share/era/data/data_series/"
    ; fils1=systemfunc("ls "+diri+"era_sur_global.nc")
    
      year1=1979
      year2=2012
      yearend=201212
      f1=addfile(fils1, "r")
      date=cd_calendar(f1->time, 1)
      timeind=ind(date.ge.197901.and.date.le.201212)
      olr_r=f1->cldc(timeind,:,:)
      olr=(/olr_r@scale_factor*olr_r+olr_r@add_offset/)
      copy_VarCoords(olr_r, olr)
      copy_VarAtts(olr_r, olr)
      olr!0="time"
      delete(olr&time)
      olr&time=date(timeind)
      printVarSummary(olr)
  
      ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
      opt=(/0,5/)
      ;;先插值lon方向的
      olr=(/linmsg_n(olr, opt, 2)/)
      ;;插值lat方向的
      opt=(/0,1/)
      olr=(/linmsg_n(olr, opt, 1)/)
      printVarSummary(olr)
     
      ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
      opt=(/0,5/)
      ;;先插值lon方向的
      olr=(/linmsg_n(olr, opt, 2)/)
      ;;插值lat方向的
      opt=(/0,1/)
      olr=(/linmsg_n(olr, opt, 1)/)
      printVarSummary(olr)
     ;;;
 


      ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
      olr_djf=month_to_season(olr, "DJF")
      ; olr_djf=month_to_annual(olr, 1)
      printVarSummary(olr_djf)
      olr_djf!0="time"
      delete(olr_djf&time)
      olr_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
      ; olr_djf=(/month_to_annual(olr, 1)/)
      olr_djf_mean=dim_avg_n_Wrap(olr_djf(:,:,:), (/1,2/))

     
    do i=0,dimsizes(olr_djf&time)-1
    olr_djf(i,:,:)=(/olr_djf(i,:,:)-olr_djf_mean(i)/)
    end do

      
    olr_rc=regCoef_n(ispan(1,dimsizes(olr_djf&time),1), olr_djf, 0, 0)
    copy_VarCoords(olr_djf(2,:,:), olr_rc)
    olr_rc=(/olr_rc*10.0/)

          ;;;区域平滑
        delete(opt)
        opt = 1
        wgt  = (/1.0, 3.0, 5.0, 3.0, 1.0/)
        wgt  = wgt/sum(wgt)             ; normalize
        olr_rc  = wgt_runave_n_Wrap(olr_rc, wgt, opt, 0)
  
        wgt2  = (/1.0, 2.0, 1.0/)
        wgt2  = wgt2/sum(wgt2)             ; normalize
        
        olr_rc   = wgt_runave_n_Wrap(olr_rc, wgt2, opt, 1)

 

     wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/ICOADS CLOUDE纬向风的趋势")
    res=True
    
    res@cnFillPalette="testcmap"
    ; res@cnFillPalette="GMT_panoply"
    res@cnFillOn=True
    res@cnLinesOn=False
    res@cnFillMode="CellFill"
    cnMissingValFillColor="white"
     res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.7
    res@cnMaxLevelValF       =   0.7
    res@cnLevelSpacingF      =  0.05
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             =0
    res@mpMaxLonF             =360
    res@mpShapeMode="FreeAspect"

    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    plot=gsn_csm_contour_map(wks, olr_rc, res)
    end 