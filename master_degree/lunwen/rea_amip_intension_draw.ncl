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
    diri="/mnt/f/share/rea_msf/"
    f1=addfile(diri+"intension_rea_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    
    inten_rea_pa=f1->intension_pa(:,timeind)
    inten_rea_id=f1->intension_id(:,timeind)
    inten_rea_alt=f1->intension_alt(:,timeind)
   
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;;read amip
    diri="/mnt/f/share/CMIP5/msf_AMIP/"
    f1=addfile(diri+"intension_amip_r1i1p1_jja.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    inten_amip_pa=f1->intension_pa(:,timeind)
    inten_amip_id=f1->intension_id(:,timeind)
    inten_amip_alt=f1->intension_alt(:,timeind)
   
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;
printVarSummary(inten_amip_id)
printVarSummary(inten_amip_pa)
corr = escorc_n(inten_amip_id, inten_amip_pa , 1 ,1)
print(corr)
exit()






models=(/"amip","ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)
;;;;


;;;anomaly
inten_rea_pa= (/runave((/dim_rmvmean_n(inten_rea_pa,1)/),3,1)/)
inten_rea_id=(/runave((/dim_rmvmean_n(inten_rea_id,1)/),3,1)/)
inten_rea_alt=(/runave((/dim_rmvmean_n(inten_rea_alt,1)/),3,1)/)

inten_amip_pa= (/runave((/dim_rmvmean_n(inten_amip_pa,1)/),3,1)/)
inten_amip_id=(/runave((/dim_rmvmean_n(inten_amip_id,1)/),3,1)/)
inten_amip_alt=(/runave((/dim_rmvmean_n(inten_amip_alt,1)/),3,1)/)

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
 

mnmx_pa=ploy(inten_rea_pa,inten_amip_pa)
mnmx_id=ploy(inten_rea_id,inten_amip_id)
mnmx_alt=ploy(inten_rea_alt,inten_amip_alt)

time=ispan(1979,2008,1)


;;;************************************************************
;画图
;;;**********************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/三个环流的强度与再分析资料的比较")
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
res@trYMaxF          =2.0
res@trYMinF          =-2.0

res@tmXTOn=False
res@tmYROn=False
res@gsnYRefLine     =0.0
res@xyMonoLineColor    = False             ; want colored lines
res@xyLineColors       = (/"Black","deeppink","Blue"/) ; colors chosen
res@xyLineThicknesses	 = (/4.,4.,3./)      ; line thicknesses
res@xyDashPatterns	 = (/0.,0.,0./)      ; make all lines solid
res@vpHeightF=0.18
res@vpWidthF=0.8
res@gsnLeftString="(a)IWCI"
mid_plot(0)= gsn_csm_xy (wks,time,mnmx_id(0:1,:),res)       ; create line plot

res@trYMaxF          =2.0
res@trYMinF          =-2.0
res@gsnLeftString="(b)PWCI"
mid_plot(1)= gsn_csm_xy (wks,time,mnmx_pa(0:1,:),res)       ; create line plot
res@trYMaxF          =1.5
res@trYMinF          =-1.5
res@gsnLeftString="(c)AWCI"
mid_plot(2)= gsn_csm_xy (wks,time,mnmx_alt(0:1,:),res)       ; create line plot

delete(res@gsnLeftString)
delete(res@xyLineColors)
res@gsnXYFillColors = "deeppink"
res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
res@gsnXYFillOpacities=0.1

amip_plot(0)= gsn_csm_xy (wks,time,mnmx_id(2:3,:),res)       ; create line plot
amip_plot(1)= gsn_csm_xy (wks,time,mnmx_pa(2:3,:),res)       ; create line plot
amip_plot(2)= gsn_csm_xy (wks,time,mnmx_alt(2:3,:),res)       ; create line plot

do i=0,2
overlay(mid_plot(i),amip_plot(i))
; overlay(mid_plot(i),cmip_plot(i))
end do

pres=True

gsn_panel(wks, mid_plot, (/3,1/), pres)
end 