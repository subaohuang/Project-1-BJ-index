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
    f1=addfile(diri+"intension_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    locationxi_rea=f1->intension_pa(:,timeind)
 
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;read amip
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"intension_amip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    locationxi_amip=f1->intension_pa(:,timeind)

    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;; read cmip
   ;;有海温的几个模式的width情况

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

 diri="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 print("dim="+dim)
 width_trend=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim, "float")
 sst_trend=asciiread("/mnt/f/share/CMIP5/sst_cha/his_sstca_trend.ascii" ,dim, "float")

 sst_poid=models(ind(sst_trend.ge.0.05))
 width_poid=models(ind(width_trend.ge.0.1))
 sst_neid=models(ind(sst_trend.le.-0.1))
 width_neid=models(ind(width_trend.le.-0.1))

 ; xi_model =sst_poid(get1Dindex(sst_poid,width_poid))
 ; dong_model =sst_neid(get1Dindex(sst_neid,width_neid))
 
 ; xi_model= width_poid
 ; dong_model=width_neid
 xi_model=(/"GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/) 
 dong_model=(/"MPI-ESM-MR","MIROC5"/)
 xi_ind=get1Dindex(models, xi_model)
 
 dong_ind=get1Dindex(models, dong_model)


diri="/mnt/f/share/CMIP5/msf_history/"
f1=addfile(diri+"intension_cmip_r1i1p1_djf.nc", "r")
date=f1->time
timeind=ind(date.ge.1979.and.date.le.2008)

locationxi_cmip=f1->intension_pa(xi_ind,timeind)
locationdong_cmip=f1->intension_pa(dong_ind,timeind)

printVarSummary(locationxi_cmip)


;;;anomaly
; locationxi_rea= (/runave((/dim_rmvmean_n(locationxi_rea,1)/),3,1)/)
; locationxi_amip= (/runave((/dim_rmvmean_n(locationxi_amip,1)/),3,1)/)
; locationxi_cmip= (/runave((/dim_rmvmean_n(locationxi_cmip,1)/),3,1)/)
; locationdong_cmip= (/runave((/dim_rmvmean_n(locationdong_cmip,1)/),3,1)/)

locationxi_rea= (/dim_rmvmean_n(locationxi_rea,1)/)
locationxi_amip= (/dim_rmvmean_n(locationxi_amip,1)/)
locationxi_cmip= (/dim_rmvmean_n(locationxi_cmip,1)/)
locationdong_cmip= (/dim_rmvmean_n(locationdong_cmip,1)/)


; ;---Convert to 1D
;   a1D      = ndtooned(locationdong_amip)
;   dsizes_a = dimsizes(locationdong_amip)

; ;---Resolve the 1D indices back to their original 3D array.
;   indices  = ind_resolve(minind(a1D),dsizes_a)
;   print(indices)

;  locationdong_amip(indices(0,0),indices(0,1))=-4
mnmx_xi=ploy(locationxi_rea,locationxi_amip,locationxi_cmip)
mnmx_dong=ploy(locationxi_rea,locationxi_amip,locationdong_cmip)


time=ispan(1979,2008,1)
printVarSummary(mnmx_xi)
printVarSummary(mnmx_dong)

;;;************************************************************
;画图
;;;**********************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_amip_cmip_不同模式")
mid_plot=new(2,graphic)
amip_plot=mid_plot
cmip_plot=mid_plot

res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False             ; don't advance frame yet

res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	         = 0.6

res@trXMaxF            = 2008
res@trXMinF	         = 1979              ; set x-axis minimum
res@trYMaxF          =2
res@trYMinF          =-2

res@tmXTOn=False
res@tmYROn=False
res@gsnYRefLine     =0.0
res@xyMonoLineColor    = False             ; want colored lines
res@xyLineColors       = (/"Black","deeppink","Blue"/) ; colors chosen
res@xyLineThicknesses	 = (/4.,4.,3./)      ; line thicknesses
res@xyDashPatterns	 = (/0.,0.,0./)      ; make all lines solid

res@gsnLeftString="(a)strong"
mid_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(0:2,:),res)       ; create line plot

res@trYMaxF          =2
res@trYMinF          =-2
res@gsnLeftString="(b)weak"
mid_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(0:2,:),res)       ; create line plot

delete(res@gsnLeftString)
delete(res@xyLineColors)
res@gsnXYFillColors = "deeppink"
res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
res@gsnXYFillOpacities=0.1

amip_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(3:4,:),res)       ; create line plot
amip_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(3:4,:),res)       ; create line plot



res@gsnXYFillColors = "blue"
res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
res@gsnXYFillOpacities=0.1
cmip_plot(0)= gsn_csm_xy (wks,time,mnmx_xi(5:6,:),res)       ; create line plot
cmip_plot(1)= gsn_csm_xy (wks,time,mnmx_dong(5:6,:),res)       ; create line plot


do i=0,1
; overlay(mid_plot(i),amip_plot(i))
overlay(mid_plot(i),cmip_plot(i))
end do

pres=True

gsn_panel(wks, mid_plot, (/2,1/), pres)
end 