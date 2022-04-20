function  ploy(locationxi_rea,locationxi_amip)
begin
mnmx      = new ((/4,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
dim1=dimsizes(locationxi_amip&time)
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_avg_n_Wrap(locationxi_amip, 0)
mnmx(2,0:dim1-1) = dim_max_n_Wrap(locationxi_amip, 0)
mnmx(3,0:dim1-1) = dim_min_n_Wrap(locationxi_amip, 0)
return mnmx
end 
;;***********************************************************************
;;主程序
;;;**********************************************************************
begin

;;;;; 读入强度
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"intension_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    intension_pa_rea = f1->intension_pa(:,timeind)
    intension_id_rea = f1->intension_id(:,timeind)
    intension_alt_rea = f1->intension_alt(:,timeind)
 
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;read amip
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"intension_amip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    intension_pa_amip = f1->intension_pa(:,timeind)
    intension_id_amip = f1->intension_id(:,timeind)
    intension_alt_amip = f1->intension_alt(:,timeind)
 
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)

    ;;;anomaly
    intension_pa_rea = (/runave((/dim_rmvmean_n(intension_pa_rea,1)/),3,0)/)
    intension_id_rea = (/runave((/dim_rmvmean_n(intension_id_rea,1)/),3,0)/)
    intension_alt_rea = (/runave((/dim_rmvmean_n(intension_alt_rea,1)/),3,0)/)

    intension_pa_amip = (/runave((/dim_rmvmean_n(intension_pa_amip,1)/),3,0)/)
    intension_id_amip = (/runave((/dim_rmvmean_n(intension_id_amip,1)/),3,0)/)
    intension_alt_amip = (/runave((/dim_rmvmean_n(intension_alt_amip,1)/),3,0)/)

    intension_id_rea = (/-1*intension_id_rea/)
    intension_id_amip = (/-1*intension_id_amip/)
;;;


;;; read  location data
    
    ;;;读入rea位置指数，计算宽度

    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"location_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    locationxi_rea=f1->locationxi(:,timeind)
    locationdong_rea=f1->locationdong(:,timeind)
    width_pa_rea = locationxi_rea 
    width_pa_rea = (/locationdong_rea - locationxi_rea/)
    width_id_rea = width_pa_rea
    width_id_rea = (/locationxi_rea - 40.0/) 

    delete(f1)
    f1=addfile(diri+"location_rea_alt_r1i1p1_djf.nc", "r")
    
    locationxi_rea=f1->locationxi(:,timeind)
    locationdong_rea=f1->locationdong(:,timeind)

    width_alt_rea = locationxi_rea 
    width_alt_rea = (/locationdong_rea - locationxi_rea/)
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
    ;;;读入amip位置指数，计算宽度

    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"location_amip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    locationxi_amip=f1->locationxi(:,timeind)
    locationdong_amip=f1->locationdong(:,timeind)

    width_pa_amip = locationxi_amip 
    width_pa_amip = (/locationdong_amip - locationxi_amip/)
    width_id_amip = width_pa_amip
    width_id_amip = (/locationxi_amip - 40.0/) 

    delete(f1)
    f1=addfile(diri+"location_alt_amip_r1i1p1_djf.nc", "r")
    
    locationxi_amip=f1->locationxi(:,timeind)
    locationdong_amip=f1->locationdong(:,timeind)

    width_alt_amip = locationxi_amip 
    width_alt_amip = (/locationdong_amip - locationxi_amip/)
    ;;;;

    ;;;width anomaly
    width_pa_rea = (/runave((/dim_rmvmean_n(width_pa_rea,1)/),3,0)/)
    width_id_rea = (/runave((/dim_rmvmean_n(width_id_rea,1)/),3,0)/)
    width_alt_rea = (/runave((/dim_rmvmean_n(width_alt_rea,1)/),3,0)/)

    width_pa_amip = (/runave((/dim_rmvmean_n(width_pa_amip,1)/),3,0)/)
    width_id_amip = (/runave((/dim_rmvmean_n(width_id_amip,1)/),3,0)/)
    width_alt_amip = (/runave((/dim_rmvmean_n(width_alt_amip,1)/),3,0)/)
    ;;;




;;;;


; ;---Convert to 1D
;   a1D      = ndtooned(locationdong_amip)
;   dsizes_a = dimsizes(locationdong_amip)

; ;---Resolve the 1D indices back to their original 3D array.
;   indices  = ind_resolve(minind(a1D),dsizes_a)
;   print(indices)

;  locationdong_amip(indices(0,0),indices(0,1))=-4
 

    mnmx_xi = ploy(intension_id_rea, intension_id_amip)
    mnmx_dong = ploy(intension_pa_rea,intension_pa_amip)
    mnmx_width = ploy(intension_alt_rea, intension_alt_amip)

    mnmx_id = ploy(width_id_rea, width_id_amip)
    mnmx_pa = ploy(width_pa_rea,width_pa_amip)
    mnmx_alt = ploy(width_alt_rea, width_alt_amip)
   



    time=ispan(1979,2008,1)
    printVarSummary(mnmx_xi)

;;;************************************************************
;画图
;;;**********************************************************
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_amip_three_circulation_intension_location")
    mid_plot=new(6,graphic)
    amip_plot=mid_plot

    res                    = True              ; plot mods desired
    res@gsnDraw            = False             ; don't draw yet
    res@gsnFrame           = False             ; don't advance frame yet

    res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
    res@vpWidthF 	         = 0.6

    res@trXMaxF            = 2008
    res@trXMinF	         = 1979              ; set x-axis minimum
    res@trYMaxF          =2.0
    res@trYMinF          =-2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@gsnYRefLine     =0.0
    res@xyMonoLineColor    = False             ; want colored lines
    res@xyLineColors       = (/"Black","deeppink","Blue"/) ; colors chosen
    res@xyLineThicknesses	 = (/4.,4.,3./)      ; line thicknesses
    res@xyDashPatterns	 = (/0.,0.,0./)      ; make all lines solid
   
    res@gsnLeftString="(a)-IWCI"
    mid_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(0:1,:),res)       ; create line plot

    res@gsnLeftString="(d)IWCW"
    res@trYMaxF          =30.0
    res@trYMinF          =-30.0
    mid_plot(1)= gsn_csm_xy (wks,time,mnmx_id(0:1,:),res)       ; create line plot

    res@trYMaxF          =2.0
    res@trYMinF          =-2.0
    res@gsnLeftString="(b)PWCI"
    mid_plot(2)= gsn_csm_xy (wks,time,mnmx_dong(0:1,:),res)       ; create line plot

    res@gsnLeftString="(e)PWCW"
    res@trYMaxF          =30.0
    res@trYMinF          =-30.0
    mid_plot(3)= gsn_csm_xy (wks,time,mnmx_pa(0:1,:),res)       ; create line plot


    res@trYMaxF          =1.6
    res@trYMinF          =-1.6
    res@gsnLeftString="(c)AWCI"
    mid_plot(4)= gsn_csm_xy (wks,time,mnmx_width(0:1,:),res)       ; create line plot
    
    res@gsnLeftString = "(f)AWCW"
    res@tiYAxisString = ""
    res@trYMaxF          =16.0
    res@trYMinF          =-16.0
    mid_plot(5)= gsn_csm_xy (wks,time,mnmx_alt(0:1,:),res)       ; create line plot


    delete(res@gsnLeftString)
    delete(res@xyLineColors)
    res@gsnXYFillColors = "deeppink"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    res@gsnXYFillOpacities=0.1

    amip_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(2:3,:),res)       ; create line plot
    amip_plot(2)= gsn_csm_xy (wks,time,mnmx_dong(2:3,:),res)       ; create line plot
    amip_plot(4)= gsn_csm_xy (wks,time,mnmx_width(2:3,:),res)       ; create line plot

    amip_plot(1)= gsn_csm_xy (wks,time,mnmx_id(2:3,:),res)       ; create line plot
    amip_plot(3)= gsn_csm_xy (wks,time,mnmx_pa(2:3,:),res)       ; create line plot
    amip_plot(5)= gsn_csm_xy (wks,time,mnmx_alt(2:3,:),res)       ; create line plot


    do i=0,5
    overlay(mid_plot(i),amip_plot(i))
    end do

    pres=True

    gsn_panel(wks, mid_plot, (/3,2/), pres)
    end 