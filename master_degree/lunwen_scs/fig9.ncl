
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
    ;;;

    print("*****************读入海温场*********************")
   ;;;海温场
    delete(date)
    delete(timeind)
    fils1=systemfunc("ls /mnt/f/share/ersst/HadIsst_sst.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    print("******************日期*********************")
    timeind=ind(date.le.201212.and.date.ge.192001) 
    sst_mre_had=lonFlip(f1->sst(timeind,{-60:60},:))
    delete(sst_mre_had&time)
    sst_mre_had&time = date(timeind)
    delete(date)
    delete(f1)
    delete(fils1)
    sst_mre_had=(/rmMonAnnCycTLL(sst_mre_had)/)
    sst_mre_had_annual=month_to_annual(sst_mre_had, 1)
    sst_mre_had_annual!0="time"
    sst_mre_had_annual&time=ispan(1920,2012,1)
    printVarSummary(sst_mre_had_annual)
    sst_annual=sst_mre_had_annual
     ;;;IPO
    sst_ipo_r=(/dim_avg_n_Wrap(sst_annual(:,{5:-5},{210:270}), (/1,2/))\
    -  (dim_avg_n_Wrap(sst_annual(:,{45:25},{150:210}), (/1,2/))+dim_avg_n_Wrap(sst_annual(:,{-25:-45},{150:210}), (/1,2/)))/2\
    /)
    ; sst_ipo_sm=runave_n_Wrap(sst_ipo, 9, 1, 0)
    ; sst_ipo_sm=(/-1*sst_ipo_sm/) 
    ; sst_ipo_sm=(/dim_standardize_n_Wrap(sst_ipo_sm, 0, 0)/)


   ;;;
    ipo_index=asciiread("/mnt/f/share/ersst/ipo_index.ascii",93,"float");;;數據是從1920年開始的 ，1979
    rc=regCoef_n(ipo_index, sst_ipo_r,0,0)
    sst_ipo=rc*sst_ipo_r+rc@yintercept
    sst_ipo_sm=runave_n_Wrap(sst_ipo(59:88), 9, 1, 0)
    sst_ipo_sm=(/-1*sst_ipo_sm/) 
  



locationxi_rea= (/runave((/dim_rmvmean_n(locationxi_rea,1)/),9,1)/)
locationxi_amip= (/runave((/dim_rmvmean_n(locationxi_amip,1)/),9,1)/)
locationdong_amip= (/runave((/dim_rmvmean_n(locationdong_amip,1)/),9,1)/)
locationxi_cmip = (/runave((/dim_rmvmean_n(locationxi_cmip,1)/),9,1)/)
locationdong_cmip= (/runave((/dim_rmvmean_n(locationdong_cmip,1)/),9,1)/)
location_amip=(/runave((/dim_rmvmean_n(location_amip,1)/),9,1)/)


sst_cha_xi=(/runave(sst_cha_xi, 9, 1)/)
sst_cha_dong=(/runave(sst_cha_dong, 9, 1)/)
sst_cha_rea=(/runave(sst_cha_rea, 9, 1)/)





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



; mnmx3      = new ((/5,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
; mnmx3(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
; mnmx3(1:,:)=sst_cha_xi


; mnmx4      = new ((/3,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
; mnmx4(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
; mnmx4(1:,:)=sst_cha_dong


mnmx3      = new ((/6,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx3(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
mnmx3(1,:)=(/sst_ipo_sm/)
mnmx3(2:,:)=sst_cha_xi


mnmx4      = new ((/4,dimsizes(locationxi_rea&time)/), typeof(locationxi_rea))
mnmx4(0,:)=dim_avg_n_Wrap(sst_cha_rea, 0)
mnmx4(1,:)=(/sst_ipo_sm/)
mnmx4(2:,:)=sst_cha_dong


 m= escorc_n(mnmx3(0,:), mnmx3(1,:), 0, 0)
 print("m="+m)

 std=dim_stddev_n(mnmx3, 1)
 print(std)
 

 

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
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/scs_fig.9")
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
res@trYMaxF          =1
res@trYMinF          =-1




res@tmXBLabelFontHeightF=0.015
res@tmXBLabelFont =  "times-roman"
res@tmYLLabelFont = "times-roman"
res@tmXBMajorLengthF=0.01        ;調整刻度的长度
res@tmXBMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用

res@tmXBMinorLengthF=0.006        ;調整刻度的长度
res@tmXBMinorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用


res@tmXTMajorLengthF=0.01        ;調整刻度的长度
res@tmXTMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用

res@tmYLMajorLengthF=0.01    
res@tmYLMajorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用    ;調整刻度的长度

res@tmYLMinorLengthF=0.006    
res@tmYLMinorOutwardLengthF=0  ;调整刻度的长度，注意2个命令要同时使用    ;調整刻度的长度

res@tmXBLabelDeltaF=-0.8
res@tmBorderThicknessF=2.0
res@tmYLLabelFontHeightF=0.015

res@tmYLMode="Explicit"
res@tmYLValues=fspan(-1.0, 1.0, 5)
res@tmYLMinorValues=fspan(-1.0, 1.0, 21)
res@tmYLLabels=(/"-1.0","-0.5","0","0.5","1.0"/)



res@gsnYRefLine     =0.0
res@xyMonoLineColor    = False             ; want colored lines
res@xyLineColors       = (/"Black","red","Blue","green","yellow","gray"/) ; colors chosen
res@xyLineThicknesses	 = (/4.,4.,4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns	 = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid

; res@gsnLeftString="(a)~F22~+PWCI"
res@tmXBLabelsOn=True
mid_plot(0)= gsn_csm_xy (wks,time,mnmx1,res)       ; create line plot

res@trYMaxF          =1
res@trYMinF          =-1
; res@gsnLeftString="(b)~F22~-PWCI"
delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
res@xyLineColors       = (/"Black","red","Blue","green"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0.,0./)      ; make all lines solid
res@tmYLLabelsOn=False
mid_plot(1)= gsn_csm_xy (wks,time,mnmx2,res)       ; create line plot

delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
; res@gsnLeftString="(c)~F5~D~F22~SST(+PWCI)"

; res@xyLineColors       = (/"Black","Blue","green","yellow","gray"/) ; colors chosen
res@xyLineColors       = (/"Black","red","Blue","green","yellow","gray"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid
res@tmYLLabelsOn=True
res@tmXBLabelsOn=False

mid_plot(2)= gsn_csm_xy (wks,time,mnmx3,res)       ; create line plot

delete(res@xyLineColors)
delete(res@xyLineThicknesses)
delete(res@xyDashPatterns)
; res@gsnLeftString="(d)~F5~D~F22~SST(-PWCI)"
res@xyLineColors       = (/"Black","red","Blue","green"/) ; colors chosen
res@xyLineThicknesses    = (/4.,4.,4.,4./)      ; line thicknesses
res@xyDashPatterns   = (/0.,0.,0.,0.,0./) 
res@tmYLLabelsOn=False

mid_plot(3)= gsn_csm_xy (wks,time,mnmx4,res)       ; create line plot



   ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"Black","red","Blue","green","yellow","gray"/) ; colors chosen
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
      textres@lgLabelOpacity = (/0/) 
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



             
         delete(lineres@lgLineColors)
         delete(textres@lgLabels)
         gres@XPosPercent = 60
         lineres@lgLineColors       = (/"red"/)
         textres@lgLabels=(/"-IPO Index"/)

          mid_plot(2) = simple_legend(wks,mid_plot(2),gres,lineres,textres)   

         delete(lineres@lgLineColors)
         delete(textres@lgLabels)
         gres@XPosPercent = 60
         lineres@lgLineColors       = (/"red"/)
         textres@lgLabels=(/"-IPO Index"/)

         mid_plot(3) = simple_legend(wks,mid_plot(3),gres,lineres,textres)   
;;;;;


;;; 为了修改为time new roman 字体，得自己修改字体

      txres = True
      txres@txFont = "times-roman"
      txres@txFontHeightF = 0.015
      txres@txJust = "CenterLeft"

      text1 = (/"ENS","AMIP(+CMIP5)","GFDL-CM3",\
                        "IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/)
      dum_text1 = gsn_add_text(wks, mid_plot(0), text1(::-1),1984.5,fspan(0.4, 0.9, 6) , txres)
      
      text2 = (/"ENS","AMIP(-CMIP5)","MPI-ESM-MR","MIROC5"/)
      dum_text2 = gsn_add_text(wks, mid_plot(1), text2(::-1),1984.5,fspan(0.6, 0.9, 4) , txres)
      
      text3 =  (/"Observation","GFDL-CM3","IPSL-CM5A-MR","INMCM4","CSIRO-Mk3-6-0"/)
      text5 = (/"-IPO Index"/)
      dum_text3 = gsn_add_text(wks, mid_plot(2), text3(::-1),1984.5,fspan(0.5, 0.9, 5) , txres)
      dum_text5 = gsn_add_text(wks, mid_plot(2), text5(::-1),1999.5,0.9 , txres)

       text4 = (/"Observation","MPI-ESM-MR","MIROC5"/)
      text6 = (/"-IPO Index"/)
      dum_text4 = gsn_add_text(wks, mid_plot(3), text4(::-1),1984.5,fspan(0.7, 0.9, 3) , txres)
      dum_text6 = gsn_add_text(wks, mid_plot(3), text6,1999.5, 0.9 , txres)






;;; add annotation a,b,c,d
txres=True
txres@txFontHeightF=0.02

lavels=(/"c","d","a","b"/)
dum=new(4,graphic)
do i=0,3
dum(i)=gsn_add_text(wks,mid_plot(i),lavels(i), 1980, 0.9, txres)
end do

  
pres=True
pres@gsnPanelXF=(/0.0416,0.5016,0.0416,0.5016/)
pres@gsnPanelDebug=True
gsn_panel(wks, (/mid_plot(2),mid_plot(3),mid_plot(0),mid_plot(1)/), (/2,2/), pres)




end 