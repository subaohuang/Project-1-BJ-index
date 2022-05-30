;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
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


function Lanczos_highpass(var)
begin
    month1 = 10*12*1.0
    ; month2 = 12.0
    nwt = 49     ;;; loose 24 months each end
    fca = 1./month1 ;;;通过的频率
    ; fcb = 1./month2
    ihp = 1   ;;;1 high pass  0:low pass         
    nsigma = 1.  ; Lanczos sigma

    wgts = filwgts_lanczos (nwt, ihp, fca, -999, nsigma)  

    var  = (/wgt_runave_n_Wrap(var, wgts,0,0)/)
 return var
 end

 function bw_bandpass(var)
 begin
    ca    = 10.0*12.0        ; band start (longer period)
    cb    = 12.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    opt   = False       ; use default options (time series of filtered
                        ; values will be returned)
    var = (/bw_bandpass_filter (var,fca,fcb,opt,0)/)  
 return var
 end 

function composite_mean(ssta_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
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
    endmon = 99712
    startyear = 1
    endyear = 997

   ; ; ;;;读入海温数据
    f1 = addfile("ssta_picontrol_highpass.nc", "r")
    ssta = f1->ssta
    ; ssta = ssta*86400.0
;;;;;;

; ; ;;;;
    f2 = addfile("uv_picontrol_highpass.nc", "r")
    u_850 = f2->u_850
    v_850 = f2->v_850
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile("piControl_monsoon_and_ensoyear_highpass.nc", "r") 
    ;;;; 03-997

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


;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"(a) El Nino","(b) La Nina","(c) Normal"/)
   str1= (/"El Nino","La Nina","Normal"/)

  plot = new(3,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "Highpass-CESM2-picontrol-ssta-evolution-during-diffmonsoon")

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
   

do k=0,2 
  print("***********start*********************"+k)

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])


    com_month_str  = composite_mean(ssta_month_12,ind_ts1)
    ssta_str_evo(0:11) = com_month_str[0]
  
    com_month_weak  = composite_mean(ssta_month_12,ind_tw1)
    ssta_weak_evo(0:11) = com_month_weak[0]
    
    ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)
   

    com_month_nor  = composite_mean(ssta_month_12,ind_tnor1_r)
   
    ssta_nor_evo(0:11) = com_month_nor[0]


    com_month_str  = composite_mean(ssta_month_12(0:5,:,:,:),ind_ts1+1)
    
    ssta_str_evo(12:) = com_month_str[0]
  

    com_month_weak  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tw1+1)
    
    ssta_weak_evo(12:) = com_month_weak[0]
    
    
    com_month_nor  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tnor1_r+1)
    
    ssta_nor_evo(12:) = com_month_nor[0]

    ;;;;写入数据
    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
    data := (/ssta_str_evo,ssta_weak_evo,ssta_nor_evo/)
    
    opt1=str1(k)+""+"diffmonsoon-SSTevolution-highpass.txt"
    ; write_matrix(data, "12f7.2",opt1)
    asciiwrite (opt1, data)

    

;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  res@tmXTOn         = False
  res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = 17

  res@gsnYRefLine = 0.0

  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 17, 1)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
  res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
   res@tmYLMode = "Explicit"

  if k.eq.2 then 
   res@trYMinF        = -1.0
   res@trYMaxF        = 1.0

   res@tmYLValues := fspan(-1.0, 1.0, 11)
   res@tmYLLabels := fspan(-1.0, 1.0, 11)
  else
   res@trYMinF        = -2.0
   res@trYMaxF        = 2.0

   res@tmYLValues := fspan(-2.0, 2.0, 11)
   res@tmYLLabels := fspan(-2.0, 2.0, 11)
  end if 

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.02


  res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1,1/)
  res@xyMarkerSizes = (/0.03,0.03,0.03/)
  res@xyMarkerColors= (/"red","blue","black"/)


  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = (/"red","blue","black"/)
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  plot(k)  = gsn_csm_xy (wks, ispan(0,17,1), (/ssta_str_evo,ssta_weak_evo,ssta_nor_evo/), res)
  ;;添加legend 
    gres = True
    if k.eq.0 then 
    gres@YPosPercent = 25.    ; expressed as %, 0->100, sets position of top border of legend 
    else
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    end if
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True

    textres@lgLabels = (/"Strong ISMR","Weak ISMR", "Normal ISMR"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
   
end do
pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,2/), pres)
end 