;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************
function spb_cal(data,data2,mxlag)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start June 
     acr = new((/3,mxlag/),float)  ;;;超前
     printVarSummary(data2)

     do i=0,2
         do j=0,mxlag-1
              acr(i,j) = escorc(data(i,0:yearnum-1),data2(i+j+1,0:yearnum-1))
         end do 
     end do
   ;;;
   
return acr
end 
;;;;
begin
  ;;;read data
  lat1 = -10
  lat2 =  30
  lon1 =  100
  lon2 = 180

  ; startmon = 198201
  ; endmon = 201912

  startmon = 197901
  endmon = 201512

  ; f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
  ; date:=cd_calendar(f1->time, 1)
  ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
  ; olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

  ; printVarSummary(olr)


  f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")

  date:=cd_calendar(f1->time, 1)
  timeind :=ind(date.ge.startmon.and.date.le.endmon)
  olr = f1->precip(timeind,{lat1:lat2},{lon1:lon2}) ;;;2.5 ×2.5 grid


  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.le.201512.and.date.ge.197901)

  sst_r =  lonFlip(f2->sst(timeind,:,:))

  sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

  printVarSummary(sst)


  fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
  f3 = addfiles(fils3, "r")

  date:=cd_calendar(f3[:]->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)
  time = f3[:]->time(timeind)
  lat = f3[1]->latitude
  lon = f3[1]->longitude


  div_r = new((/dimsizes(time), dimsizes(lat),dimsizes(lon)/),float)


  do i = 0,dimsizes(timeind)/12-1

  div_r(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d(:,{200},:,:))

  end do
 
  ; div = dim_avg_n_Wrap(div_r(:,:,{lat1-1:lat2+1},{lon1-1:lon2+1}),1)
  div = div_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})


  div = (/div*10^5/)
  div = (/-1*div/)

  ; div = (/rmMonAnnCycTLL(div)/)


;;;;chazhi 
   print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    xo = olr&lon
    yo = olr&lat



    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo


    div_chazhi = linint2_Wrap(div&longitude, div&latitude, div, True, xo, yo, 0)

    div_chazhi!1= "lat"
    div_chazhi!2= "lon"

    div_chazhi&lat = yo
    div_chazhi&lon = xo
;;;;
   delete(sst)
   delete(div)
   
;;;;
    printVarSummary(sst_chazhi)
    olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;;

;;;;anomaly 

    olr = (/rmMonAnnCycTLL(olr)/)
    sst_chazhi = (/rmMonAnnCycTLL(sst_chazhi)/)
    div_chazhi = (/rmMonAnnCycTLL(div_chazhi)/)
;;;

;;;;Monthly bin

   ; sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   ; olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   dim = dimsizes(sst_chazhi)
   
  sst_chazhi_june = sst_chazhi(5::12,:,:)
  sst_chazhi_july = sst_chazhi(6::12,:,:)
  sst_chazhi_aug = sst_chazhi(7::12,:,:)

  sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
  sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


  olr_june = olr(5::12,:,:)
  olr_july = olr(6::12,:,:)
  olr_aug = olr(7::12,:,:)

  olr_jj =  array_append_record(olr_june,olr_july,0)
  olr_jja = array_append_record(olr_jj, olr_aug,0)


  ; div_june = div_chazhi(2::12,:,:)
  ; div_july = div_chazhi(3::12,:,:)
  ; div_aug = div_chazhi(4::12,:,:)

  ; div_jj =  array_append_record(div_june,div_july,0)
  ; div_jja = array_append_record(div_jj , div_aug,0)

   ;  sst_one = ndtooned(sst_chazhi_jja(:,{0:15},{125:160}))
   ;  olr_one = ndtooned(olr_jja(:,{0:15},{125:160}))

    ; sst_one = ndtooned(sst_chazhi(:,{0:20},{125:160}))
;;;;

;;;;计算区域平均

   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi_jja(:,{0:20},{125:160}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr_jja(:,{0:20},{125:160}), (/1,2/))


   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi_jja(:,{5:25},{125:150}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr_jja(:,{5:25},{125:10}), (/1,2/))


   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi(:,{0:20},{125:160}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr(:,{0:20},{125:160}), (/1,2/))

    sst_wnp = dim_avg_n_Wrap(sst_chazhi(:,{0:20},{125:160}), (/1,2/))

    olr_wnp = dim_avg_n_Wrap(olr(:,{0:20},{125:160}), (/1,2/))

  ; sst_wnp = dim_avg_n_Wrap(sst_chazhi(:,{5:25},{120:150}), (/1,2/))

  ; olr_wnp = dim_avg_n_Wrap(olr(:,{5:25},{120:150}), (/1,2/))

   ; sst_wnp = dim_avg_n_Wrap(sst_chazhi(:,{5:25},{120:150}), (/1,2/))

   ; olr_wnp = dim_avg_n_Wrap(olr(:,{5:25},{120:150}), (/1,2/))
; ;;


;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(sst_wnp)/12 -1
  sst_wnp_12 = new((/12,yearnum/),float)
  olr_wnp_12 = sst_wnp_12

  month = 11-5   ;;5,6,7
  do i=0,11   
    if i.le.month then 
     sst_wnp_12(i,:) = sst_wnp(i+11-month:dimsizes(sst_wnp)-1-12:12)
     olr_wnp_12(i,:) = olr_wnp(i+11-month:dimsizes(olr_wnp)-1-12:12)
     else
     sst_wnp_12(i,:) = sst_wnp(i+11-month:dimsizes(sst_wnp)-1:12)
     olr_wnp_12(i,:) = sst_wnp(i+11-month:dimsizes(olr_wnp)-1:12)
    end if
  end do 

  ;;;时间都是1979-6 -2015-05


;;;计算超前滞后回归
    
    x = sst_wnp
    y = olr_wnp

    mxlag    = 8
    
  
    ccr = new ( (/3, 2*mxlag+1/), float)    
    

    x_Lead_y = spb_cal(sst_wnp_12,olr_wnp_12,mxlag)     ;;SST 超前 OLR 

    y_Lead_x = spb_cal(olr_wnp_12,sst_wnp_12,mxlag)     ;;SST 滞后 OLR 


   
    ccr(:,0:mxlag-1) = x_Lead_y(:,0:mxlag-1:-1)     ; "positive lag"

    ccr(:,mxlag) = escorc_n(sst_wnp_12(0:2,:),olr_wnp_12(0:2,:),1,1)
    
    ccr(:, mxlag+1:) = y_Lead_x(:,0:mxlag-1) 


    ; acr_lead = new((/6/),float)  ;;;超前
    ; acr_lag = new((/6/),float)  ;;;超前
    
    ;  do j=0,5
    ;       acr_lead(j) = escorc(sst_wnp_12(0,0:yearnum-1),olr_wnp_12(j+1,0:yearnum-1))
    ;  end do 
    
     acr_test = escorc(sst_wnp_12(2, 0:yearnum-1),olr_wnp_12(8,0:yearnum-1))

     print(acr_test)
     exit()
     print(ccr(:,mxlag))

    

    wks=gsn_open_wks("pdf", "WNP-JJA-PRE&hadisst超前滞后相关")
    ; wks=gsn_open_wks("pdf", "SCS-JJA-OLR&hadisst超前滞后相关")

    res=True
    res@xyLineColors = (/"black","red","blue"/)
    res@xyDashPatterns = (/0,0,0/)
    res@gsnYRefLine = 0.0
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    res@tmXBMinorOn = False
    res@tmXTMinorOn = False
    res@tmXBValues = ispan(-mxlag,mxlag,1)
    res@trXMinF = -mxlag
    res@trXMaxF = mxlag
    res@trYMaxF = 1.0
    res@trYMinF = -1.0

    plot=gsn_csm_xy(wks, ispan(-mxlag,mxlag,1), ccr, res)

    ; plot=gsn_csm_xy(wks, fspan(-6,6,13), dim_avg_n_Wrap(ccr, 0) , res)
  

 end 



