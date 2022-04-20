
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
    
    location_amip=f1->intension_pa(:,timeind)
  
    delete(timeind)
    delete(date)
    delete(f1)
    delete(diri)
;;;; read cmip
   ;;有海温的几个模式的width情况

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

 xi_model=(/"GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/) 
 dong_model=(/"MPI-ESM-MR","MIROC5"/)

 xi_ind=get1Dindex(models, xi_model)
 dong_ind=get1Dindex(models, dong_model)

locationxi_amip=location_amip(xi_ind,:)
locationdong_amip=location_amip(dong_ind,:)


 diri="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 print("dim="+dim)
 width_trend=asciiread(diri+"cmip_msf_inten_trend.ascii" ,dim, "float")

    ;;；读入海温差的时间序列
    diri="/mnt/f/share/CMIP5/sst_cha/"
    f1=addfile(diri+"sstca_cmip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)
    sst_cha=f1->sst_cha(:,timeind)
   


    diri="/mnt/f/share/CMIP5/msf_history/"
    f1=addfile(diri+"intension_cmip_r1i1p1_djf.nc", "r")
    date=f1->time
    timeind=ind(date.ge.1979.and.date.le.2008)

    locationxi_cmip=f1->intension_pa(xi_ind,timeind)
    locationdong_cmip=f1->intension_pa(dong_ind,timeind)

    sst_cha_xi=sst_cha(xi_ind,:)
    sst_cha_dong=sst_cha(dong_ind,:)

    printVarSummary(locationxi_cmip)
    delete(f1)
    delete(diri)

;;;读入观测的海温数据

    ;;；读入海温差的时间序列
    diri="/mnt/f/share/ersst/"
    f1=addfile(diri+"sstca_rea_r1i1p1_djf.nc", "r")
    sst_cha_rea=f1->sst_cha(:,timeind)
   

;;;anomaly
; locationxi_rea= (/runave((/dim_rmvmean_n(locationxi_rea,1)/),3,1)/)
; locationxi_amip= (/runave((/dim_rmvmean_n(locationxi_amip,1)/),3,1)/)
; locationxi_cmip= (/runave((/dim_rmvmean_n(locationxi_cmip,1)/),3,1)/)
; locationdong_cmip= (/runave((/dim_rmvmean_n(locationdong_cmip,1)/),3,1)/)
;;;
 

locationxi_rea= (/dim_rmvmean_n(locationxi_rea,1)/)
locationxi_amip= (/dim_rmvmean_n(locationxi_amip,1)/)
locationdong_amip= (/dim_rmvmean_n(locationdong_amip,1)/)
locationxi_cmip = (/dim_rmvmean_n(locationxi_cmip,1)/)
locationdong_cmip= (/dim_rmvmean_n(locationdong_cmip,1)/)
location_amip=(/dim_rmvmean_n(location_amip,1)/)



mnmx1      = new ((/6,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx1(0,:)=dim_avg_n_Wrap(locationxi_rea, 0)
mnmx1(1,:)=dim_avg_n_Wrap(locationxi_amip, 0)
mnmx1(2:,:)=locationxi_cmip


mnmx2      = new ((/4,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx2(0,:)=dim_avg_n_Wrap(locationxi_rea, 0)
mnmx2(1,:)=dim_avg_n_Wrap(locationdong_amip, 0)
mnmx2(2:,:)=locationdong_cmip
;;;SST 

sst_cha_rea= (/dim_rmvmean_n(sst_cha_rea,1)/)

sst_cha_xi = (/dim_rmvmean_n(sst_cha_xi,1)/)
sst_cha_dong= (/dim_rmvmean_n(sst_cha_dong,1)/)



mnmx3      = new ((/5,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx3(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
mnmx3(1:,:)=sst_cha_xi


mnmx4      = new ((/3,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx4(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
mnmx4(1:,:)=sst_cha_dong


 m= escorc_n(mnmx1(0,:), mnmx3(0,:), 0, 0)
 print("m="+m)
 ; exit()

; ;---Convert to 1D
;   a1D      = ndtooned(locationdong_amip)
;   dsizes_a = dimsizes(locationdong_amip)

; ;---Resolve the 1D indices back to their original 3D array.
;   indices  = ind_resolve(minind(a1D),dsizes_a)
;   print(indices)

;  locationdong_amip(indices(0,0),indices(0,1))=-4


time=ispan(1979,2008,1)


;;;************************************************************
;画图
;;;**********************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/rea_amip_cmip_不同模式")
mid_plot=new(4,graphic)
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
res@xyLineColors       = (/"Black","red","Blue","green","yellow","gray"/) ; colors chosen
res@xyLineThicknesses	 = (/4.,4.,4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns	 = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid

res@gsnLeftString="(a)~F22~+PWCI"
mid_plot(0)= gsn_csm_xy (wks,time,mnmx1,res)       ; create line plot

res@trYMaxF          =2
res@trYMinF          =-2
res@gsnLeftString="(b)~F22~-PWCI"
delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
res@xyLineColors       = (/"Black","red","Blue","green"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0.,0./)      ; make all lines solid
mid_plot(1)= gsn_csm_xy (wks,time,mnmx2,res)       ; create line plot

delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
res@gsnLeftString="(c)~F5~D~F22~SST(+PWCI)"

res@xyLineColors       = (/"Black","Blue","green","yellow","gray"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0.,0./)      ; make all lines solid


mid_plot(2)= gsn_csm_xy (wks,time,mnmx3,res)       ; create line plot

delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
res@gsnLeftString="(d)~F5~D~F22~SST(-PWCI)"
res@xyLineColors       = (/"Black","Blue","green"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0./) 

mid_plot(3)= gsn_csm_xy (wks,time,mnmx4,res)       ; create line plot



   ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"Black","red","Blue","green","yellow","gray"/) ; colors chosen
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"ENS","AMIP(+CMIP5)","GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/)
   
    mid_plot(0) = simple_legend(wks,mid_plot(0),gres,lineres,textres)   
    delete(lineres@lgLineColors)
    delete(textres@lgLabels)
    lineres@lgLineColors       = (/"Black","red","Blue","green"/) ; colors chosen
     textres@lgLabels=(/"ENS","AMIP(-CMIP5)","MPI-ESM-MR","MIROC5"/)
     mid_plot(1) = simple_legend(wks,mid_plot(1),gres,lineres,textres)   

     delete(lineres@lgLineColors)
     delete(textres@lgLabels)
     lineres@lgLineColors       = (/"Black","Blue","green","yellow","gray"/) ; colors chosen
     textres@lgLabels=(/"Observation","GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/)
      mid_plot(2) = simple_legend(wks,mid_plot(2),gres,lineres,textres)   


        delete(lineres@lgLineColors)
        delete(textres@lgLabels)
        lineres@lgLineColors       = (/"Black","Blue","green"/) ; colors chosen
        textres@lgLabels=(/"Observation","MPI-ESM-MR","MIROC5"/)
        mid_plot(3) = simple_legend(wks,mid_plot(3),gres,lineres,textres)   


  
pres=True
gsn_panel(wks, mid_plot, (/2,2/), pres)




end 