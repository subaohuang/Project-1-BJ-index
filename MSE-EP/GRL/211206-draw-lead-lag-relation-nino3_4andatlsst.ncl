function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12
    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 


function corr_lead_lag(Nino3_4_djf,ssta_12,lat1,lat2,lon1,lon2)
begin
    
    atl = dim_avg_n_Wrap(ssta_12(:,:,{lat1:lat2},{lon1:lon2}),(/2,3/))
    dim = dimsizes(atl)

    atl_3runmean = runave_n_Wrap(atl(:,0:dim(1)-2),3, 0, 0)
    atl_3runmean(0:dim(0)-3,:) = atl_3runmean(1:dim(0)-2,:) ;;;;10 number JFM
    atl_3runmean(dim(0)-2,:) =  (atl(10,0:dim(1)-2) + atl(11,0:dim(1)-2) + atl(0,1:dim(1)-1))/3                ;;;add NDJ
    atl_3runmean(dim(0)-1,:) =  (atl(11,0:dim(1)-2) + atl(0,1:dim(1)-1) + atl(1,1:dim(1)-1))/3                ;;;add NDJ
 
    corr_lead = new(20,float)
    dims = dimsizes(atl_3runmean)
    corr_lead(0:11) = escorc_n(Nino3_4_djf(1:dim(1)-2), atl_3runmean(0:11,0:dims(1)-2), 0, 1) ;;JFM DJF corr 
    corr_lead(12:) = escorc_n(Nino3_4_djf(1:dim(1)-2), atl_3runmean(0:7,1:dims(1)-1), 0, 1)

return corr_lead
end 


function part_corr(var1,var2,var3)
begin 
  ;;;remove var3 influence 
  r12 = escorc_n(var1, var2, 0, 0)
  ; copy_VarCoords(var2(0,:,:), r12)  ;;;这是二维
  r13 = escorc_n(var1, var3, 0, 0)  ;;;这个是1dim
  r23 = escorc_n(var2, var3, 0, 0)   ;;;这个是2dim
  r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
  rc_mam = r12
  rc_mam = (/r12_3/)
return rc_mam
end 

function band_pass(data)
begin
  ihp   = 2                            ; band_pass
  sigma = 1.0                           ; Lanczos sigma
  nWgt  = 9                           ; will lose nwgt-1 data must odd number 
  fcb   = 1./2.                        ; 2 years
  fca   = 1./9.                         ;7年
  wgts  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
  msf_500_band=wgt_runave_Wrap(data, wgts, 0)
  return msf_500_band
  end 
begin
    ;;;time span
    startmon =198001
    endmon = 202012

    running_time = 15

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;; 
;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon-100)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta = sst
    ssta=(/rmMonAnnCycTLL(sst)/)


    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    copy_VarCoords(sst, ssta)

    ssta_12 = split_12(ssta) ;;; month,year,lat,lon
  ;;
  

    Nino3_4 = dim_avg_n_Wrap(ssta(:,{-5:5},{190:240}),(/1,2/))
    Nino3_4_djf = month_to_season(Nino3_4, "DJF") ;;;1980 JF , 1980D 1981JF
    
    corr_atl = corr_lead_lag(Nino3_4_djf,ssta_12,10,30,280,330)
    corr_ind = corr_lead_lag(Nino3_4_djf,ssta_12,-30,30,40,100)
    corr_nino34 = corr_lead_lag(Nino3_4_djf,ssta_12,-5,5,190,240)
   
  
    plot = new(1,graphic)
;;;;

    ; rc_mam_u = where(prob_mam_u.le.0.05.or.prob_mam_v.le.0.05, rc_mam_u, rc_mam_u@_FillValue)
    ; rc_mam_v = where(prob_mam_u.le.0.05.or.prob_mam_v.le.0.05, rc_mam_v, rc_mam_u@_FillValue)
 ; wks=gsn_open_wks("pdf", "lead-lag-corr-EAESM-nino3.4-ATLrun")
    wks=gsn_open_wks("pdf", "lead-lag-corr-ATL-Nino34")


    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@xyLineColors     = (/"black","red","blue"/)          ; change line color
    res@xyDashPatterns = (/0,0,0/)

    res@xyMonoMarker = True
    res@xyMarker = 1

    res@xyLineThicknessF=3.0
  
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF = 1.0
    res@trYMinF = -1.0
    res@tmYLValues    = ispan(-10,10,2)*0.1
    res@tmYLLabels    = ispan(-10,10,2)*0.1
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
   
   
    res@trXMaxF = 19
    res@trXMinF = 0
    res@tmXBValues    = ispan(2, 19, 3)
    res@tmXBMinorValues = ispan(0, 19, 1)
    res@tmXBLabels    = (/"MAM(-1)","JJA(-1)","SON(-1)","D(0)JF(+1)","MAM(+1)","JJA(+1)","SON(+1)"/)

    res@gsnXRefLine = 17
    res@gsnYRefLine = 0


    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.008  ;调整刻度的长度，注意2个命令要同时使用
    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8
    res@tmBorderThicknessF=2.0

    res@tmYLLabelFontHeightF=0.015
    ; res@vpXF=0.15
    ; res@vpYF=0.9
    res@vpHeightF=0.4
    res@vpWidthF=0.6
    
    ; res@tiMainString=" Spring Double ITCZ Bias and DJF Enino Intensity"
    res@tiMainString=" EAESM  and Nino3 lead-lag-corr"
    ; res@tiMainFontHeightF = 0.015
    ; res@tiYAxisString="~F22~ SSTA std change[k]"
    ; res@tiYAxisFontHeightF=0.012
    
    ; res@tiXAxisString="~F22~ dsst change [k]  "

    ; res@tiXAxisFontHeightF=0.012
    
    plot=gsn_csm_xy(wks, ispan(0,19,1), (/corr_nino34,corr_atl,corr_ind/), res)


    ;;legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    lineres@lgLineColors = res@xyLineColors ; line colors
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"Nino3_4","North Atlantic","Indian"/)

    plot = simple_legend(wks,plot,gres,lineres,textres)   

    draw(plot)
    frame(wks)
end 