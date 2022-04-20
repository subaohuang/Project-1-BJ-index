;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 

function lead_corr(fils1,fils2)
begin

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

   
    time1 = f1->time
    ; time1 = (/time1 -1/)

    date1=cd_calendar(time1, 1)

    timeind1=ind(date1.ge.197901.and.date1.le.201612)
   
    ;;;;U,V
    u_CT_200  = short2flt(f1->u(timeind1,{200},:,:))
    u_CT_850  =  short2flt(f1->u(timeind1,{850},:,:))

    u_CT_850 = (/rmMonAnnCycTLL(u_CT_850)/)
    u_CT_200 = (/rmMonAnnCycTLL(u_CT_200)/)

    u_CT_200 =  (/dtrend_msg_n(ispan(1,dimsizes(u_CT_200&time),1),u_CT_200,False,True,0)/)
    u_CT_850 =  (/dtrend_msg_n(ispan(1,dimsizes(u_CT_850&time),1),u_CT_850,False,True,0)/)



;;;;
   
    time2 = f2->time
    ; time2 = (/time2 -1/)

    date2=cd_calendar(time2, 1)

    timeind2=ind(date2.ge.197901.and.date2.le.201612)
   
  
    sst_ctrl = lonFlip(f2->sst(timeind2,:,:))
    

    sst_ctrl = (/rmMonAnnCycTLL(sst_ctrl)/)


    sst_ctrl =  (/dtrend_msg_n(ispan(1,dimsizes(sst_ctrl&time),1),sst_ctrl,False,True,0)/)


; ;;;; Nino3
;     Nino3_ctrl = dim_avg_n_Wrap(sst_ctrl(:,{-5:5},{210:270}), (/1,2/))
;     dim = dimsizes(Nino3_ctrl)
;     print(dim)
; ;;；

; ;;;; Nino3.4 170w-120w
;     Nino3_ctrl = dim_avg_n_Wrap(sst_ctrl(:,{-5:5},{190:240}), (/1,2/))
;     dim = dimsizes(Nino3_ctrl)
;     print(dim)
; ;;; 
    


;;;; TASST -30-30, 290-360
    Nino3_ctrl = dim_avg_n_Wrap(sst_ctrl(:,{-5:5},{290:360}), (/1,2/))
    dim = dimsizes(Nino3_ctrl)
    print(dim)
;;; 



; ;;;;;DMI指数
;     DMI_ctrl = u_CT_850(:,1,1)
  
;     DMI_ctrl  =  (/dim_avg_n_Wrap(u_CT_850(:,{5:20},{40:110}), (/1,2/)) - \
;                   dim_avg_n_Wrap(u_CT_200(:,{5:20},{40:110}), (/1,2/))/)
;     DMI_ctrl_JJA = month_to_season(DMI_ctrl(24:dim-13), "JJA")  ;;;jja 少一年，从25601开始
; ;;;;


;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_CT_850(:,1,1)
  
  
    DMI_ctrl  =  (/dim_avg_n_Wrap(u_CT_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_CT_850(:,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl(24:dim-13), "JJA")  ;;;jja 少一年，从25601开始
;;;;


;;;;Nino3 season 
  Nino3_ctrl_JJA_ne1 =  month_to_season(Nino3_ctrl(12:dim-25), "JJA") ;;;25506 ,07,08
 
  Nino3_ctrl_SON_ne1 =  month_to_season(Nino3_ctrl(12:dim-25), "SON") 

  Nino3_ctrl_DJF = month_to_season(Nino3_ctrl, "DJF") ;; 25401-25402:0 25412-25501-25502:1,25512-25601-25602:2 

  Nino3_ctrl_DJF_0 =  Nino3_ctrl_DJF(2:dim/12-2)      ;25512 25601


  Nino3_ctrl_MAM_0 =  month_to_season(Nino3_ctrl(24:dim-13), "MAM") ;;;25603 ,04,05
 
  Nino3_ctrl_JJA_0 =  month_to_season(Nino3_ctrl(24:dim-13), "JJA") 
 
  
  Nino3_ctrl_SON_0 =  month_to_season(Nino3_ctrl(24:dim-13), "SON") 
 
  Nino3_ctrl_DJF_1 =  Nino3_ctrl_DJF(3:dim/12-1)      ;25612 25701
 
  Nino3_ctrl_MAM_1 =  month_to_season(Nino3_ctrl(36:dim-1), "MAM") ;;;25703 ,04,05
  
  Nino3_ctrl_JJA_1 =  month_to_season(Nino3_ctrl(36:dim-1), "JJA") ;;;25703 ,04,05
 
;;;
 
;;;;相关系数


   JJAne1_ctrl = escorc(Nino3_ctrl_JJA_ne1,DMI_ctrl_JJA)
 

   SONne1_ctrl = escorc(Nino3_ctrl_SON_ne1,DMI_ctrl_JJA)


   DJF0_ctrl = escorc(Nino3_ctrl_DJF_0,DMI_ctrl_JJA)
  

   MAM0_ctrl = escorc(Nino3_ctrl_MAM_0,DMI_ctrl_JJA)


   JJA0_ctrl = escorc(Nino3_ctrl_JJA_0,DMI_ctrl_JJA)

   SON0_ctrl = escorc(Nino3_ctrl_SON_0,DMI_ctrl_JJA)


   DJF1_ctrl = escorc(Nino3_ctrl_DJF_1,DMI_ctrl_JJA)
  

   MAM1_ctrl = escorc(Nino3_ctrl_MAM_1,DMI_ctrl_JJA)
  
   JJA1_ctrl = escorc(Nino3_ctrl_JJA_1,DMI_ctrl_JJA)

   corr_ctrl = (/JJAne1_ctrl,SONne1_ctrl, DJF0_ctrl, MAM0_ctrl, JJA0_ctrl, SON0_ctrl, DJF1_ctrl, MAM1_ctrl,JJA1_ctrl/)
   
   return corr_ctrl
end 

begin
    ;;;cam4 输出结果
    fils1 = "/home/yangsong3/data-observation/ERA-interim/monthly/uwnd.mon.mean.nc"
    ;;;pop 输出结果
    fils2 = "/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc"
    
     
    
    corr_ATLyear_ctrl = lead_corr(fils1,fils2)


   ; wks=gsn_open_wks("pdf", "lead-lag-corr-EAESM-nino3.4-ERAI")
   wks=gsn_open_wks("pdf", "lead-lag-corr-EAESM-TASST-ERAI")

    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@xyLineColors     = (/"black"/)          ; change line color
    res@xyDashPatterns = (/0/)

    res@xyMonoMarker = True
    res@xyMarker = 1

    res@xyLineThicknessF=3.0
  
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF = 0.8
    res@trYMinF = -0.8
    res@tmYLValues    = ispan(-8,8,2)*0.1
    res@tmYLLabels    = ispan(-8,8,2)*0.1
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
   
   
    res@trXMaxF =  8
    res@trXMinF = 0
    res@tmXBValues    = ispan(0, 8, 1)
    res@tmXBLabels    = (/"JJA(-1)","SON(-1)","DJF(0)","MAM(0)","JJA(0)","SON(0)","DJF(1)","MAM(1)","JJA(1)"/)

    res@gsnXRefLine = 4
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
    res@tiMainString=" EASM  and Nino3.4 lead-lag-corr"
    ; res@tiMainFontHeightF = 0.015
    ; res@tiYAxisString="~F22~ SSTA std change[k]"
    ; res@tiYAxisFontHeightF=0.012
    
    ; res@tiXAxisString="~F22~ dsst change [k]  "

    ; res@tiXAxisFontHeightF=0.012
    
    plot=gsn_csm_xy(wks, ispan(0,8,1), (/corr_ATLyear_ctrl/), res)


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
    textres@lgLabels = (/"ERAI"/)

    plot = simple_legend(wks,plot,gres,lineres,textres)   


    draw(plot)
    frame(wks)
end 