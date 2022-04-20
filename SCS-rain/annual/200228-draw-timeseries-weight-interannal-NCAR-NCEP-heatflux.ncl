;;********************************************************
;;计算SST倾向以及SST值与净长波辐射，短波辐射，感热，潜热以及平流和对流,
;;其中热通量的方向都是向上为正，向下为负。。
;;;的回归场，热通量数据是NCEP-NCAR数据，海洋数据是SODA2.2.4,SST是哈哈得来SST数据
;;由于数据之间的精度不同，需要进行插值预处理,SODA精度是0.5×0.5,HADISST是1×1
;;NCEP-NCAR的精度为1.87×1.9
;;s
;;                                         lnshuheng 2020.02.18
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
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160
 

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_annual = month_to_annual(ssta(:,{lat1:lat2},{lon1:lon2}), 1)
    ssta_annual =(/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)

    ;;;sst 倾向值
     dsst = center_finite_diff_n(ssta, 1, False, 0, 0) 
     copy_VarCoords(ssta, dsst)
     dsst_annual = month_to_annual(dsst(:,{lat1:lat2},{lon1:lon2}), 1)
     dsst_annual =(/dtrend_msg_n(ispan(1,yearnum,1),dsst_annual,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    ;;;dsst calculate by annual mean method

    dsst_new = center_finite_diff_n(ssta_annual, 1, False, 0, 0) 
    
    copy_VarCoords(ssta_annual, dsst_new)
    printVarSummary(dsst_new)

    dsst_new2 = dsst_new
    dsst_new2(0:33,:,:) = ssta_annual(1:34,:,:)-ssta_annual(0:33,:,:)


    ;;;;





;;;

;;get region mean as index     
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   ssta_index = dim_avg_n_Wrap(ssta_annual(:,{0:15},{110:155}), (/1,2/))
   dsst_index = dim_avg_n_Wrap(dsst_annual(:,{0:15},{110:155}), (/1,2/))
   dsst_new_index = dim_avg_n_Wrap(dsst_new(:,{0:15},{110:155}), (/1,2/))
   dsst_new2_index = dim_avg_n_Wrap(dsst_new2(:,{0:15},{110:155}), (/1,2/))

   ;;;normalized 
   
   ssta_index = (/dim_standardize_n(ssta_index, 1, 0)/)
   dsst_index = (/dim_standardize_n(dsst_index, 1, 0)/)
   dsst_new_index = (/dim_standardize_n(dsst_new_index, 1, 0)/)
   dsst_new2_index = (/dim_standardize_n(dsst_new2_index, 1, 0)/)
;;;



;;;calculate dsst_








;;;read NCEP-NCAR heat flux data
    
 
    fils2 = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
    f2 = addfile(fils2, "r")

    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)
  
    
    lhtfl = f2->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f2->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f2->nlwrs(timeind,{lat1:lat2},{lon1:lon2})
    
   
    lhtfl  = (/rmMonAnnCycTLL(lhtfl)/)
    shtfl  = (/rmMonAnnCycTLL(shtfl)/)
    nswrs  = (/rmMonAnnCycTLL(nswrs)/)
    nlwrs  = (/rmMonAnnCycTLL(nlwrs)/)

    ;;;转化成方向向下

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)

   lhtfl_a = dim_avg_n_Wrap(lhtfl(:,{0:15},{110:155}), (/1,2/))
   shtfl_a = dim_avg_n_Wrap(shtfl(:,{0:15},{110:155}), (/1,2/))
   nswrs_a = dim_avg_n_Wrap(nswrs(:,{0:15},{110:155}), (/1,2/))
   nlwrs_a = dim_avg_n_Wrap(nlwrs(:,{0:15},{110:155}), (/1,2/))
  
   

    lhtfl_annual = month_to_annual(lhtfl_a, 1)
    shtfl_annual = month_to_annual(shtfl_a, 1)
    nswrs_annual = month_to_annual(nswrs_a, 1)
    nlwrs_annual = month_to_annual(nlwrs_a, 1)
    
    total_annual = lhtfl_annual
    total_annual = lhtfl_annual + shtfl_annual + nswrs_annual +nlwrs_annual
    
    
    lhtfl_annual =(/dtrend_msg_n(ispan(1,yearnum,1),lhtfl_annual,False,True,0)/)
    shtfl_annual =(/dtrend_msg_n(ispan(1,yearnum,1),shtfl_annual,False,True,0)/)
    nswrs_annual =(/dtrend_msg_n(ispan(1,yearnum,1),nswrs_annual,False,True,0)/)
    nlwrs_annual =(/dtrend_msg_n(ispan(1,yearnum,1),nlwrs_annual,False,True,0)/)
    total_annual =(/dtrend_msg_n(ispan(1,yearnum,1),total_annual,False,True,0)/)


    lhtfl_annual = (/dim_standardize_n(lhtfl_annual, 1, 0)/)
    shtfl_annual = (/dim_standardize_n(shtfl_annual, 1, 0)/)
    nswrs_annual = (/dim_standardize_n(nswrs_annual, 1, 0)/)
    nlwrs_annual = (/dim_standardize_n(nlwrs_annual, 1, 0)/)
    total_annual = (/dim_standardize_n(total_annual, 1, 0)/)
   ;;;

;;;;weight total flux using montly data to get annual data 

   weight1 = ispan(1,12,1)*1.0
   weight1 = (2*weight1-1)/24
   
   weight2 = ispan(1,12,1)*1.0
   weight2 = (25-2*weight2)/24
   
   total_month = lhtfl_a
   total_month = lhtfl_a + shtfl_a + nswrs_a + nlwrs_a
   
   total_annual_test = total_annual
   do i =0,34,1
      k := mod(i+1,2)
       if k.eq.0 then 
         total_annual_test(i) = sum(total_month(i*12:(i+1)*12-1)*weight2)
       else 
         total_annual_test(i) = sum(total_month(i*12:(i+1)*12-1)*weight1)
        end if
   end do

   printVarSummary(total_annual_test)


  total_annual_test =(/dtrend_msg_n(ispan(1,yearnum,1),total_annual_test,False,True,0)/)
  total_annual_test = (/dim_standardize_n(total_annual_test, 1, 0)/)
;;;;;   






           


;;;;



;;;;
  

print("***********画时间序列***************")

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "dSST-test&DSST以及总热通量各项的时间序列-test")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2013
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1979,2013,5)
  resxy@tmXBLabels = ispan(1979,2013,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2013,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknesses =  (/ 5,3,2,2,2,2,2/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"black","black","red","red"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,1,0,1/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines"/)
  
  time = ispan(1979, 2013, 1)
  plot = gsn_csm_xy(wks,time,(/dsst_new_index, dsst_new2_index,total_annual_test,total_annual/), resxy)

  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = resxy@xyLineColors
  lineres@lgLineThicknesses = (/ 5,5,5,5/)
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns


  textres = True
  textres@lgLabels = (/"dSST-test","DSST","total_annual_weight","total_annual"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)

  end
