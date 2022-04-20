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

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
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

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta = sst
    ssta=(/rmMonAnnCycTLL(sst)/)


    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    copy_VarCoords(sst, ssta)
   ;;
;;;
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_ncep1_50-20.nc", "r")

    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    startmon = 198001
    endmon = 202012

    startyear = 1980
    endyear = 2020

    timeind :=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)
   

    ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)
   

    ; MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_nep = MSE_col
    MSE_nep = (/MSE_col/10^6/)


;;;;;;split to 12 month  is a array size (12,yearnum)
  
     MSE_nep_12 = split_12(MSE_nep)
     ssta_12 = split_12(ssta)
   
;;;;;

   
;;;等到年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual!0 = "time"
    MSE_nep_annual&time = ispan(startyear,endyear,1)

    ssta_annual = month_to_annual(ssta, 1)
    ssta_annual!0 = "time"
    ssta_annual&time = ispan(startyear,endyear,1)
    
    ssta_annual = (/dim_avg_n(ssta_12(0:4,:,:,:),0)/)
    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)

    ssta_mam = ssta_annual
;;;;;

;;;
 
    MSE_nep_mam = dim_avg_n_Wrap(MSE_nep_annual(:,{0:20},{250:280}), (/1,2/))



    ;;;局地海温
    ; ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

    ;  ;;nino1+2:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{270:280}), (/1,2/))
    ; ;;;
    
    ;  ;;nino3:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{150:270}), (/1,2/))
    ; ;;;
     
     ;;;nino4 
      nino3_4 = dim_avg_n_Wrap(ssta_mam(:,{-5:5},{190:240}), (/1,2/))
      ; atlantic = dim_avg_n_Wrap(ssta_mam(:,{10:30},{280:330}), (/1,2/))
      atlantic = dim_avg_n_Wrap(ssta_mam(:,{0:20},{40:100}), (/1,2/))  ;;;indian
    ;;;

    ;;;;去掉局地海温的影响,偏相关性
    ;;; 1 for MSE, 2 for global sst , 3 for nino3_4 ,4 for atlantic SST
        r12no3 = part_corr(MSE_nep_mam,ssta_mam, nino3_4) ;;;12.3  2dim
        ; rc_mam_rmaltantic = part_corr(MSE_nep_mam,ssta_mam,nino3_4) ;;;12.4 
        r12no4 = part_corr(MSE_nep_mam,ssta_mam, atlantic) ;;; 2dim

        r13no4 = part_corr(MSE_nep_mam,nino3_4, atlantic)  ;;; 1dim

        r23no4 = part_corr(ssta_mam,nino3_4,atlantic)    ;;; 1dim
       
        ; rc_mam = (/r12no3/)

        r13no4_23no4 = r13no4*r23no4

        r13no4_23no4_2dim  = r12no4
        r13no4_23no4_2dim  = r13no4_23no4
        rc_mam = (/r12no4 - r13no4_23no4_2dim/)/(sqrt((1- r13no4^2)*(1- r23no4^2)))

        copy_VarCoords(ssta_mam(0,:,:), rc_mam)


    ;;;;


    ;;;;有效自由度
        ; r1_auto = esacr_n(MSE_nep_mam,1,0)
        ; r2_auto = esacr_n(ssta_mam,1,0)
        ; one = r2_auto(0,:,:)
        ; one =1
        ; N = dimsizes(ssta_mam&time)*(one- r1_auto(1)*r2_auto(1,:,:))/(one+r1_auto(1)*r2_auto(1,:,:))
        ; N_int = toint(N)
        N = dimsizes(ssta_mam(:,1,1))
        ; prob_mam = rtest(rc_mam, N-1, 0)  ;;;1个变量去掉，加1
         prob_mam = rtest(rc_mam, N-2, 0)  ;;;1个变量去掉，加1
        copy_VarCoords(ssta_mam(0,:,:), prob_mam)
    ;;;;;
    

    ; ;;;normal corr

    ;   rc_mam = escorc_n(MSE_nep_mam, ssta_mam, 0, 0)

    ;   copy_VarCoords(ssta_mam(1,:,:), rc_mam)

    ;   prob_mam = rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
    ;   copy_VarCoords(ssta_mam(0,:,:), prob_mam)
    ; ;;;;


    plot = new(4,graphic)
;;;;

    
wks   = gsn_open_wks("pdf", "MSE与SST的相关-JFMAM-remove_Nino3.4_rmindian")
; wks   = gsn_open_wks("pdf", "MSE与SST的相关-JFMAM-remove_Nino3.4_rmatlantic")
; wks   = gsn_open_wks("pdf", "MSE与SST的相关-MAM-去掉nino4的影响")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2
  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  res@lbLabelStride = 1


  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"

  ; res@gsnLeftString="(PCOR MAM MSE, MAM SST)" 
  res@gsnLeftString="Partical corr (NEP MSE & SST , remove Nino3.4 & IND SST)" 
  ; res@gsnLeftString="Partical corr (NEP MSE & SST , remove Nino3.4 & NTA SST)" 
  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam ,res)
  ; res@gsnLeftString="(b)JJA MSE trend " 
  ; plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  ; res@gsnLeftString="(c)SON MSE trend " 
  ; plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  ; res@gsnLeftString="(d)DJF MSE trend " 
  ; plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 


   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 250
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
  do i =0,3 
   plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 


 
   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(4,graphic)
     i=0
    iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)
  
  
   
     overlay(plot(i), iplot(i))
  ; end do

    ;;;  
   
; pres = True
; pres@gsnPanelLabelBar  = True
; pres@pmLabelBarOrthogonalPosF = -0.01
draw(plot(0))
frame(wks)
; gsn_panel(wks, plot, (/2,2/), pres)

end


