function  ploy(locationxi_rea,locationxi_amip,locationxi_cmip)
begin
mnmx      = new ((/7,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
dim1=dimsizes(locationxi_amip&time)
dim2=dimsizes(locationxi_cmip&time)
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_avg_n_Wrap(locationxi_amip, 0)
mnmx(2,0:dim2-1) = dim_avg_n_Wrap(locationxi_cmip, 0)
mnmx(3,0:dim1-1) = dim_max_n_Wrap(locationxi_amip, 0)
mnmx(4,0:dim1-1) = dim_min_n_Wrap(locationxi_amip, 0)
mnmx(5,0:dim2-1) = dim_max_n_Wrap(locationxi_cmip, 0)
mnmx(6,0:dim2-1) = dim_min_n_Wrap(locationxi_cmip, 0)
return mnmx
end 
;;***********************************************************************
;;主程序
;;;**********************************************************************
begin
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"location_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    locationxi_rea=f1->locationxi(:,timeind)
    locationdong_rea=f1->locationdong(:,timeind)
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;read amip
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"location_amip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    locationxi_amip=f1->locationxi(:,timeind)
    locationdong_amip=f1->locationdong(:,timeind)
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;; read cmip
    xi_model=(/"FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-ESM2M","INMCM4"/)

;;有海温的几个模式的width情况
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FIO-ESM","FGOALS-s2","FGOALS-g2","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

   xi_ind=get1Dindex(models, xi_model)
    diri="/mnt/f/share/CMIP5/msf_history/"
    f1=addfile(diri+"location_cmip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2005)
    
    locationxi_cmip=f1->locationxi(xi_ind,timeind)
    locationdong_cmip=f1->locationdong(xi_ind,timeind)
   printVarSummary(locationxi_cmip)

;;;;
width_rea=new((/dimsizes(locationxi_rea&model),dimsizes(locationxi_rea&time)/),typeof(locationxi_rea))

do i=0,dimsizes(locationxi_rea&model)-1
width_rea(i,:)= locationdong_rea(i,:) - locationxi_rea(i,:)
end do

width_amip=new((/dimsizes(locationxi_amip&model),dimsizes(locationxi_amip&time)/),typeof(locationxi_rea))
do i=0,dimsizes(locationxi_amip&model)-1
width_amip(i,:)=locationdong_amip(i,:) - locationxi_amip(i,:)
end do

width_cmip=new((/dimsizes(locationxi_cmip&model),dimsizes(locationxi_cmip&time)/),typeof(locationxi_rea))
do i=0,dimsizes(locationxi_cmip&model)-1
width_cmip(i,:)=locationdong_cmip(i,:) - locationxi_cmip(i,:)
end do
copy_VarCoords(locationxi_rea, width_rea)
copy_VarCoords(locationxi_amip, width_amip)
copy_VarCoords(locationxi_cmip, width_cmip)

;;;anomaly
locationxi_rea= (/runave((/dim_rmvmean_n(locationxi_rea,1)/),3,1)/)
locationdong_rea=(/runave((/dim_rmvmean_n(locationdong_rea,1)/),3,1)/)
width_rea=(/runave((/dim_rmvmean_n(width_rea,1)/),3,1)/)

locationxi_amip= (/runave((/dim_rmvmean_n(locationxi_amip,1)/),3,1)/)
locationdong_amip=(/runave((/dim_rmvmean_n(locationdong_amip,1)/),3,1)/)
width_amip=(/runave((/dim_rmvmean_n(width_amip,1)/),3,1)/)


; locationxi_cmip= (/runave((/dim_rmvmean_n(locationxi_cmip,1)/),3,1)/)
; locationdong_cmip=(/runave((/dim_rmvmean_n(locationdong_cmip,1)/),3,1)/)
; width_cmip=(/runave((/dim_rmvmean_n(width_cmip,1)/),3,1)/)


; ;---Convert to 1D
;   a1D      = ndtooned(locationdong_amip)
;   dsizes_a = dimsizes(locationdong_amip)

; ;---Resolve the 1D indices back to their original 3D array.
;   indices  = ind_resolve(minind(a1D),dsizes_a)
;   print(indices)

;  locationdong_amip(indices(0,0),indices(0,1))=-4
 



mnmx_xi=ploy(locationxi_rea,locationxi_amip,locationxi_cmip)
mnmx_dong=ploy(locationdong_rea,locationdong_amip,locationdong_cmip)
mnmx_width=ploy(width_rea,width_amip,width_cmip)
time=ispan(1979,2008,1)
printVarSummary(mnmx_xi)

;;;************************************************************
;画图
;;;**********************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_amip_cmip_lon")
mid_plot=new(3,graphic)
amip_plot=mid_plot
cmip_plot=mid_plot

res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False             ; don't advance frame yet

res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	         = 0.7

res@trXMaxF            = 2008
res@trXMinF	         = 1979              ; set x-axis minimum
res@trYMaxF          =30
res@trYMinF          =-30

res@tmXTOn=False
res@tmYROn=False
res@gsnYRefLine     =0.0
res@xyMonoLineColor    = False             ; want colored lines
res@xyLineColors       = (/"Black","deeppink","Blue"/) ; colors chosen
res@xyLineThicknesses	 = (/4.,4.,3./)      ; line thicknesses
res@xyDashPatterns	 = (/0.,0.,0./)      ; make all lines solid
res@vpHeightF=0.18
res@vpWidthF=0.8
res@gsnLeftString="(a)Western Edge"
mid_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(0:2,:),res)       ; create line plot

res@trYMaxF          =12
res@trYMinF          =-12
res@gsnLeftString="(b)Eastern Edge"
mid_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(0:2,:),res)       ; create line plot
res@trYMaxF          =30
res@trYMinF          =-30
res@gsnLeftString="(c)Width"
mid_plot(2)= gsn_csm_xy (wks,time,mnmx_width(0:2,:),res)       ; create line plot

delete(res@gsnLeftString)
delete(res@xyLineColors)
res@gsnXYFillColors = "deeppink"
res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
res@gsnXYFillOpacities=0.1

amip_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(3:4,:),res)       ; create line plot
amip_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(3:4,:),res)       ; create line plot
amip_plot(2)= gsn_csm_xy (wks,time,mnmx_width(3:4,:),res)       ; create line plot


res@gsnXYFillColors = "blue"
res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
res@gsnXYFillOpacities=0.1
cmip_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(5:6,:),res)       ; create line plot
cmip_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(5:6,:),res)       ; create line plot
cmip_plot(2)= gsn_csm_xy (wks,time,mnmx_width(5:6,:),res)       ; create line plot

do i=0,2
overlay(mid_plot(i),amip_plot(i))
; overlay(mid_plot(i),cmip_plot(i))
end do

pres=True

gsn_panel(wks, mid_plot, (/3,1/), pres)
end 