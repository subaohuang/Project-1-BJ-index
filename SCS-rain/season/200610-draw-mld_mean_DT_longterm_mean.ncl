;;;画ERA 的0-12小时forecast月平均的通量数据的趋势
;;
;;
;;                                         lnshuheng 2020.03.10
;**********************************************************
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
;;;;

function  wgt_annual(total_month,season)
begin

    if season.eq."MAM"
       h = 3
    else if season.eq."JJA"
        h=6
       else if season.eq."SON"
          h = 9
      else 
          h = 0
       end if 
    end if
    end if
   print(h)
    weight1 = ispan(h,h+2,1)*1.0
    weight1 = (weight1-h+0.5)/3.0

    weight2 = ispan(h,h+2,1)*1.0
    weight2 = (h-weight2+2.5)/3.0
    
    dim = dimsizes(total_month)
    printVarSummary(total_month)



    weight1_3d = conform_dims(dimsizes(total_month(1:3,:,:)), weight1, 0)
    weight2_3d = conform_dims(dimsizes(total_month(1:3,:,:)), weight2, 0)


    total_annual_test = new((/dim(0)/12,dim(1),dim(2)/),"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    h = h-1 ;;;考虑到数组index要减1  

    if season.ne."DJF" then
        do i =0,32,1
          total_annual_test(i+1,:,:) =  dim_sum_n((total_month(h+i*12:h+2+i*12,:,:)*weight1_3d),0) + dim_sum_n(total_month(h+(i+1)*12:h+2+(i+1)*12,:,:)*weight2_3d,0) + \
                                  dim_sum_n(total_month(h+3+i*12:(i+1)*12-1,:,:),0) +  dim_sum_n(total_month((i+1)*12:h-1+(i+1)*12,:,:),0)
        end do    
    end if 

    if season.eq."DJF" then
      ;;;12月的情况为
      DT12 = total_annual_test
      DT1 = DT12
      DT2 = DT12
      do i =0,31,1
          DT12(i+1,:,:) = 0.5*total_month(i*12+11,:,:) + 0.5*total_month((i+1)*12+11,:,:) + dim_sum_n(total_month((i+1)*12:(i+1)*12+10,:,:),0)
          DT1(i+1,:,:) = 0.5*total_month((i+1)*12,:,:) + 0.5*total_month((i+2)*12,:,:) + dim_sum_n(total_month((i+1)*12+1:(i+1)*12+11,:,:),0)
          DT2(i+1,:,:) = 0.5*total_month((i+1)*12+1,:,:) + 0.5*total_month((i+2)*12+1,:,:) + dim_sum_n(total_month((i+1)*12+2:(i+1)*12+11,:,:),0) + total_month((i+2)*12,:,:)
          total_annual_test(i+1,:,:) = (DT12(i+1,:,:)+DT1(i+1,:,:)+DT2(i+1,:,:))/3.0
      end do    
     total_annual_test(33,:,:) =  default_fillvalue("float")
    end if
    

return total_annual_test
end    
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   ; lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 160
;;;;



;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ; ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season)
  

    ;;;sst 倾向值
     dsst = center_finite_diff_n(ssta, 1, False, 0, 0) 
     copy_VarCoords(ssta, dsst)
     ; dsst_annual = month_to_season(dsst(:,{lat1:lat2},{lon1:lon2}), season)
     

    ;;;
    delete(date)
    delete(timeind)
;;;


   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/SODA3.4.2_mld_density_criterion_0.8.nc", "r") 
    


    indyear = (2013-1980+1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r = f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    mlt = f6->mld(0:indyear,{lat1:lat2},{lon1:lon2})
    
    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)
 
   copy_VarCoords(u_r, u) 
   copy_VarCoords(v_r, v)
   copy_VarCoords(wt_r, wt)
   copy_VarCoords(T_r, T)


    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(T_r)


 
    ;;;calculate Horizontal and Vertical advction

     
    ;;;;垂直积分平均，在混合层内进行
     z = u&depth 
    ; dz_4d = conform_dims(dimsizes(u), dz, (/1,2,3/))
      dz_4d = new(dimsizes(u),float)
      do i = 0,dimsizes(mlt&time)-1
        dz_4d(i,:,:,:) = dpres_plevel(z, mlt(i,:,:), 0, 0)
      end do 
    ;;;1混合层平均的
     T_mean = T(:,1,:,:)
     u_mean = u(:,1,:,:)
     v_mean = v(:,1,:,:)

     T_mean = (/dim_sum_n_Wrap(T*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)/)
     u_mean = (/dim_sum_n_Wrap(u*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)/)
     v_mean = (/dim_sum_n_Wrap(v*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)/)
    ;;;

    ;;;温度的梯度

    T_mam_annual = month_to_season(T_mean, "MAM")
    T_jja_annual = month_to_season(T_mean, "JJA")
    T_son_annual = month_to_season(T_mean, "SON")
    T_djf_annual = month_to_season(T_mean, "DJF")
    T_annual = month_to_annual(T_mean, 1)

    
    dT_mam  = T_mam_annual
    dT_jja  = T_mam_annual
    dT_son  = T_mam_annual
    dT_djf  = T_mam_annual

    dT_mam = (/center_finite_diff_n(T_mam_annual, 1, False, 0, 0)/)
    dT_jja = (/center_finite_diff_n(T_jja_annual, 1, False, 0, 0)/)
    dT_son = (/center_finite_diff_n(T_son_annual, 1, False, 0, 0) /)
    dT_djf = (/center_finite_diff_n(T_djf_annual, 1, False, 0, 0)/) 
    dT_annual = (/center_finite_diff_n(T_annual, 1, False, 0, 0)/) 


    dt_mam_mean = dim_avg_n_Wrap(dT_mam, 0)
 

    dt_mam_mean&lat@units = "degree_north"
    dt_mam_mean&lon@units = "degree_east"

    dt_jja_mean = dt_mam_mean
    dt_son_mean = dt_mam_mean
    dt_djf_mean = dt_mam_mean
    dt_annual_mean = dt_mam_mean

    dt_jja_mean = (/dim_avg_n_Wrap(dT_jja, 0)/)
    dt_son_mean = (/dim_avg_n_Wrap(dT_son, 0)/)
    dt_djf_mean = (/dim_avg_n_Wrap(dT_djf, 0)/)


    dt_annual_mean = (/dim_avg_n_Wrap(dT_annual, 0)/)


    dt_mam_mean = (/dt_mam_mean - dt_annual_mean/)
    dt_jja_mean = (/dt_jja_mean - dt_annual_mean/)
    dt_son_mean = (/dt_son_mean - dt_annual_mean/)
    dt_djf_mean = (/dt_djf_mean - dt_annual_mean/)




;;;;

print("*********开始画图********")
;;**************************************************************
;画图
;;***************************************************************
  ; wks   = gsn_open_wks("pdf", "所有通量的线性趋势-ERAI-newthod-mld_0.8-"+season+"wgt-mld")
  wks   = gsn_open_wks("pdf", "混合层温度DT-long-term-mean-diff-annual")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  

  ; res@cnFillPalette = "NCV_jaisnd"
   res@cnFillPalette = "BlueWhiteOrangeRed"

  res@gsnSpreadColorStart = 130
  res@gsnSpreadColorEnd = 253

  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn= False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             =lon1
  res@mpMaxLonF             = lon2

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   0.00
  res@cnMaxLevelValF       =   0.06
  res@cnLevelSpacingF      =  0.005

    res@cnMinLevelValF       =   -0.01
  res@cnMaxLevelValF       =   0.01
  res@cnLevelSpacingF      =  0.001


  res@cnInfoLabelOn = False

   

  plot = new(4,graphic)

  res@gsnLeftString = "MAM DT long-term mean "
   
  plot(0) = gsn_csm_contour_map(wks, dt_mam_mean, res)
  res@gsnLeftString = "JJA DT long-term mean "
  plot(1) = gsn_csm_contour_map(wks,dt_jja_mean,res)
  res@gsnLeftString = "SON DT long-term mean "
  plot(2) = gsn_csm_contour_map(wks, dt_son_mean, res)
  res@gsnLeftString = "DJF DT long-term mean "
  plot(3) = gsn_csm_contour_map(wks, dt_djf_mean, res)


pres = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)

end