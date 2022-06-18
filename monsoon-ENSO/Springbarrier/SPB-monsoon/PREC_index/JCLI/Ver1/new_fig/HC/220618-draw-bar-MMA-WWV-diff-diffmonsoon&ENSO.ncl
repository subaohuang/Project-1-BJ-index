;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_bw(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
function composite_mean(ssta_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 120
   lon2 = 280
   ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
    ; u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
    ; v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,{lat1:lat2},{lon1:lon2}),(/1,2,3/))

return [/ssta_month_str/]
end 


;;*************************************************
;; main program 
;;*************************************************************
begin
  

    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

    ;  startmon = 20101
    ; endmon = 120012
    ; startyear = 201
    ; endyear = 1200
   
   diri_in = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"
    ; year_out = addfile("../60-100_piControl_monsoon_and_ensoyear.nc", "r")
    year_out = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/piControl_monsoon_and_ensoyear_highpass_bw.nc", "r")

    nino_year = year_out->nino_year 
    lanina_year = year_out->lanina_year 
    normal_year = year_out->normal_year 

    monstryear = year_out->monstryear 
    monwekyear = year_out->monwekyear
    noryear = year_out->noryear
 
    wet_warm_year =  year_out->wet_warm_year
    wet_cold_year = year_out->wet_cold_year
    wet_norm_year = year_out->wet_norm_year

    dry_warm_year = year_out->dry_warm_year
    dry_cold_year = year_out->dry_cold_year
    dry_norm_year = year_out->dry_norm_year

    nor_warm_year = year_out->nor_warm_year
    nor_cold_year = year_out->nor_cold_year
    nor_norm_year = year_out->nor_norm_year
    print(num( wet_warm_year))
    print( num( wet_cold_year))  
    print( num(wet_norm_year))
    print( num(dry_warm_year))
    print( num(dry_cold_year))
    print(num(dry_norm_year))
    print(num(nor_warm_year))
    print(num(nor_cold_year))
    print(num(nor_norm_year))



   
    f3 = addfile(diri_in + "/WWV/highpass_bw_thetao_0-300mean_picontrol.nc", "r")
    time = f3->time 
    date:=cd_calendar(time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    ssta = f3->thetao(timeind,:,:)
    

;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"(a) El Nino","(b) La Nina","(c) Normal"/)
   str1= (/"El Nino","La Nina","Normal"/)

  ntilde  = "n~H-13V2F35~D~FV-2H3~"
;;;;;

do k=0,2 
  print("***********start*********************"+k)

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

;;;split to 12 month 
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 -1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  ssta_str_evo = new((/18/),float)
  ssta_weak_evo = new((/18/),float)
  ssta_nor_evo = new((/18/),float)

 do i=0,11   
     if i.le.9 then
       ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1-12:12,:,:)   ;;;3月开始，到次年2月
     else
        ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1:12,:,:)
     end if
 end do 

    com_month_str  = composite_mean(ssta_month_12,ind_ts1)
    ssta_str_evo(0:11) = com_month_str[0]
  
    com_month_weak  = composite_mean(ssta_month_12,ind_tw1)
    ssta_weak_evo(0:11) = com_month_weak[0]
    
    ; ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)
   

    com_month_nor  = composite_mean(ssta_month_12,ind_tnor1)
   
    ssta_nor_evo(0:11) = com_month_nor[0]


    com_month_str  = composite_mean(ssta_month_12(0:5,:,:,:),ind_ts1+1)
    
    ssta_str_evo(12:) = com_month_str[0]
  

    com_month_weak  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tw1+1)
    
    ssta_weak_evo(12:) = com_month_weak[0]
    
    
    com_month_nor  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tnor1+1)
    
    ssta_nor_evo(12:) = com_month_nor[0]
;;;;;;
    wwv_str = avg(ssta_str_evo(0:2))
    wwv_weak = avg(ssta_weak_evo(0:2))
    wwv_nor = avg(ssta_nor_evo(0:2))
;;;;

;;;bar 图
    wks = gsn_open_wks("pdf",""+str1(k)+"pIcontrol-BAR-MAM-WWV")          ; send graphics to PNG file
     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.3
     ; sres@vpXF = .15
     sres@trXMinF = 0.5
     sres@trXMaxF = 4.5

     ; sres@trYMinF = -1.2
     ; sres@trYMaxF = 1.2

     if k.ne.2 then 
     sres@trYMinF = -0.5
     sres@trYMaxF = 0.5
     else 
     sres@trYMinF = -0.3
     sres@trYMaxF = 0.3
     end if 

     sres@gsnDraw = False
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.3           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = ispan(1, 3, 1)
     sres@tmXBLabels = (/"Str","Normal","Weak"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "[~S~o~N~C~N~/month]" 
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
    

     str = ispan(1,2,1)
     ; color = (/"black","gray","yellow","green","pink","deeppink3","red","blue"/)
     color = (/"red"/)
     
    
     sres@gsnXYBarChartColors = color

     sres@gsnLeftStringFontHeightF = 0.02
     ; sres@gsnLeftString = "AMJJASOND mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino3 "
     sres@gsnLeftString = "(a) WWV for El Ni"+ntilde+"o"
     ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(:,i),sres)   ; draw each time series

     plot1 = gsn_csm_xy(wks,ispan(1,3,1),(/wwv_str,wwv_nor,wwv_weak/),sres)   ; draw each time series     
     color = (/"black"/)
     ; sres@gsnXYBarChartColors = color
     ; plot2 = gsn_csm_xy(wks,ispan(1,3,1),hb_mean_mamjja(1,:),sres)   ; draw each time series     

     ; sres@gsnXYBarChartColors = (/"blue"/)
     ; plot3 = gsn_csm_xy(wks,ispan(1,3,1)+0.3,hb_mean_mamjja(2,:),sres)   ; draw each time series     
    
draw(plot1)
; draw(plot2)
; draw(plot3)
frame(wks)
end do 
end 