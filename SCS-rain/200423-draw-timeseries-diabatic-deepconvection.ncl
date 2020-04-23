;;********************************************************
;;画CFSR的diabaticdeepconvection的时间序列图

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
    endmon = 200812

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region

    f2=addfile("/home/ys17-19/lsh/Project/SCS-rain/diabatic_and_deepcon_heating_CFSR_ANN_SCS.nc", "r")
  
    deepcon_q_time = f2->deepcon_q_time(:,{100000:20000})
    diab_q_time = f2->diab_q_time(:,{100000:20000})
    lev = f2->lev({100000:20000})
 
    ; printMinMax(diab_q_time, 1)
    deepcon_q_time = deepcon_q_time*86400.0
    diab_q_time = diab_q_time*86400.0
  

    ;;;垂直平均
    dp=dpres_plevel(lev,100000,min(lev),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(deepcon_q_time),dp,1)          ; confirm dp to caculate sum(var*dp)

    deepcon_q_wgt = (/dim_sum_n(dp_cfm*deepcon_q_time,1)/dim_sum(dp)/)               ; sum(var*dp)
    diab_q_wgt = (/dim_sum_n(dp_cfm*diab_q_time,1)/dim_sum(dp)/)               ; sum(var*dp)
;;;

  
print("***********画时间序列***************")

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "CFSR-diabatic-and-deepconvection的时间序列")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF =  endyear
  resxy@trXMinF = startyear
  resxy@tmXBValues = ispan(startyear,endyear,5)
  resxy@tmXBLabels = ispan(startyear,endyear,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(startyear,endyear,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknesses =  (/ 5,3,2,2,2,2,2/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"gray","black","blue","red","green","orange","black"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,1,1,1,1,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines","lines","MarkLines"/)
  
  time = ispan(startyear,endyear,1)
  plot = gsn_csm_xy(wks,time,(/diab_q_wgt,deepcon_q_wgt/), resxy)



  ; ;;;添加两条趋势线
  ; res = True
  ; res@gsLineColor = "red"
  
  ;  dum1 = gsn_add_polyline(wks, plot, ispan(1979,1997,1), line1, res)
  ;  dum2 = gsn_add_polyline(wks, plot, ispan(1998,2013,1), line2, res)


  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = (/"gray","black","blue","red","green","orange","black"/)  
  lineres@lgLineThicknesses = (/ 5,3,2,2,2,2,2/)                     ; line thicknesses
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = (/0,0,1,1,1,1,0/)


  textres = True
  textres@lgLabels = (/"GPCP"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)

  end
