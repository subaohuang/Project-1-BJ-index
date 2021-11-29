;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
function spb_cal(data,data2)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start June 
     acr = new((/12,13/),float)
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(data(i,0:yearnum-1),data(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(data(i,0:yearnum-1),data2(i+j-12,0:yearnum-1))
            end if 
         end do 
     end do
   ;;;
 ;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"
   maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   SPB = dim_sum_n_Wrap(acr_gra_max, 0)
   acr@maind = maind
   acr@SPB = SPB
return acr
end 
function split_12(nino34)
begin
    yearnum = dimsizes(nino34)/12
    nino34_12 = new((/12,yearnum/),float)

    ; do i=0,11   
    ; if i.le.6 then 
    ; nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
    ; else
    ; nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
    ; end if
    ; end do 
  do i=0,11   
    nino34_12(i,:) = nino34(i::12)
  end do 
return nino34_12
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin

 ;;;time span
    startmon =25601
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;


;;;读入控制实验组的海温数据
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl = f1->TEMP(timeind,0,:,:)
    ssta_ct = (/rmMonAnnCycTLL(sst_ctrl)/)
    copy_VarCoords(sst_ctrl, ssta_ct)

    ;;detrend 
    ssta_ct = (/dtrend_msg_n(ispan(1,dimsizes(ssta_ct&time),1),ssta_ct,False,True,0)/)
    ;;;
    nino34_ct =  dim_avg_n_Wrap(ssta_ct(:,{5:-5},{190:240}), (/1,2/))
 
;;;读入B22000 WNPCLM海温数据
    
    model_name = (/"HEAT_B_BP"/)
    ; model_name =(/"lsh_B2000_WNP_heating_Idl"/)
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_IDLthree/lsh_B2000_WNP_heating_IDLthree.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+model_name+"/"+model_name+".TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_year_80_20/lsh_B2000_alt_north_year_80_20.TEMP.h1.0251-0280.nc"
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_wnp = f2->TEMP(timeind,0,:,:)

    ; sst_wnp_r = f2->TEMP(timeind,0,:,:)
    ; sst_wnp = PopLatLon(sst_wnp_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    
    ssta_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
    copy_VarCoords(sst_wnp, ssta_wnp)

    ;;detrend 
    ssta_wnp = (/dtrend_msg_n(ispan(1,dimsizes(ssta_wnp&time),1),ssta_wnp,False,True,0)/)
    nino34_wnp =  dim_avg_n_Wrap(ssta_wnp(:,{5:-5},{190:240}), (/1,2/))
;;;;;
   
    nino34_ct_12 = split_12(nino34_ct)
    nino34_wnp_12 = split_12(nino34_wnp)

;;;;
    std_ct = dim_stddev_n(nino34_ct_12, 1)
    std_wnp = dim_stddev_n(nino34_wnp_12, 1)

 a = (/1,4,8/)
 std_1 = dim_stddev_n(a,0)
 std_2 = dim_stddev_n(a*2, 0)
 std_3 = dim_stddev_n(a+10, 0)

 print(std_1)
 print(std_2)
 print(std_3)
 exit()

plot = new(1,graphic)
wks = gsn_open_wks("pdf", "STD-NINO3-SEASONAL")
; wks1 = gsn_open_wks("pdf", "ISM-after-root_mean_square-newmonsoon.pdf")

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
  res@trXMaxF        = 11
  res@gsnYRefLine = 0.0

  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 11, 1)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
  res@tmXBLabels = "" + ispan(1, 12, 1)
  res@tmYLMode = "Explicit"


  res@trYMinF        = 0
  res@trYMaxF        = 2.0

  res@tmYLValues := fspan(0, 2., 11)
  res@tmYLLabels := fspan(0, 2., 11)


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

  plot  = gsn_csm_xy (wks, ispan(0,11,1), (/std_wnp,std_ct/), res)


  ;;添加legend 
    gres = True

    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 

                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"red","blue"/)  ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True

    textres@lgLabels = (/"Heating","ctrl"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    midplot = simple_legend(wks,plot,gres,lineres,textres)   

    draw(plot)
    frame(wks)
   end 
