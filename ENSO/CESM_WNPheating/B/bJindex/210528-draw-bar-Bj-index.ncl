
function trend(data)
begin
  rc1 =regCoef_n(ispan(1, dimsizes(data&time),1), data,0,0)
  rc1 = rc1*dimsizes(data&time)
  return rc1
end 

begin

    yearstart = 199101
    yearend = 201912

    lat1 = 0
    lat2 = 20
    lon1 = 240
    lon2 = 280

;;;;read data

    ; diri1 = "/home/ys17-19/lsh/Project/ENSO/CESM_WNPheating/B/term_BJ_CESM_ctrl.nc"
    diri1 = "/home/ys17-19/lsh/Project/ENSO/CESM_WNPheating/B/bJindex/test_term_BJ_CESM_B2000_f19g16_CP_CTRL.nc"
    ; diri1 = "/home/ys17-19/lsh/Project/ENSO/CESM_WNPheating/B/test_term_BJ_CESM_lsh_B2000_WNP_heating_Idl.nc"

    f1 = addfile(diri1, "r") 
    term_BJ_ct = f1->term_BJ(0,:)

    ; model_name = (/"lsh_B2000_WNP_heating_double"/)
    ; model_name = (/"HEAT_B_BP"/)
    model_name = (/"lsh_B2000_WNP_heating_IDLthree"/)
    ; model_name = (/"lsh_B2000_WNP_heating_Idl"/)
    ; model_name = (/"MCheating_double"/)
    diri2 = "/home/ys17-19/lsh/Project/ENSO/CESM_WNPheating/B/bJindex/test_term_BJ_CESM_"+model_name+".nc"


    f2 = addfile(diri2, "r") 
    term_BJ_wnpheating = f2->term_BJ(0,:)


    number = new((/2,6/),float)
    number(0,0:4) = term_BJ_ct 
    number(1,0:4) = term_BJ_wnpheating

    number(0,5) = sum(term_BJ_ct)
    number(1,5) = sum(term_BJ_wnpheating)

 
;==============================================================
; Mission 5. bar plot
; ==============================================================

    wks = gsn_open_wks("pdf",""+model_name+"-BJ-indexbar-test")          ; send graphics to PNG file

    sres = True
    sres@vpWidthF = 0.7
    sres@vpHeightF = 0.5
    sres@vpXF = .15


    sres@trXMinF = 0.5
    sres@trXMaxF = 6.5
    sres@trYMinF = -2.0
    sres@trYMaxF = 2.0

    sres@gsnDraw = True
    sres@gsnFrame = False
    sres@gsnXYBarChart = True
    sres@gsnXYBarChartBarWidth = 0.2           ; change bar widths
    sres@tmXBMode          = "Explicit"         ; explicit labels
    sres@tmXBValues        = (/1,2,3,4,5,6/)
    sres@tmXBLabels = (/"MA","TD","ZA","EK","TH","SUM"/)
    ; sres@tmYLPrecision = 1
    ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
    ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
    ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01

    sres@tmXBLabelFontHeightF = 0.0205
    sres@tmXTLabelFontHeightF = 0.0205
    sres@tmYLLabelFontHeightF = 0.0225
    sres@tiMainFontHeightF = 0.025
    sres@gsnRightString = ""
    sres@tiYAxisString = "yr~S~-1"
    sres@tiYAxisFontHeightF  =0.02
    sres@gsnYRefLine = 0.      ; set the reference line equal to 0


    str = ispan(1,6,1)

    color = (/"slateblue4","blue","plum","green","red","black"/)


    plot = new(2,graphic)
    ; do i = 0,1
    ; sres@gsnXYBarChartColors = color
    ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i,4.75+0.15*i,5.75+0.15*i/),number(i,:),sres)   ; draw each time series
    ; end do

    sres@vpHeightF = 0.2
    sres@vpWidthF = 0.7

    sres@vpXF = 0.2
    sres@vpYF = 0.8

    sres@gsnXYBarChartColors2 = color
    sres@gsnLeftString = "(a) Control"
    plot(0) = gsn_csm_xy(wks,ispan(1, 6, 1),number(0,:),sres)   ; draw each time series
    sres@gsnLeftString = "(b)"+ model_name +"heating minus Control"
    sres@vpYF = 0.45
    sres@trYMinF = -0.4
    sres@trYMaxF = 0.4
    sres@tmYLLabels = ispan(-4,4,2)*0.1 +""
    sres@tmYLValues = ispan(-4,4,2)*0.1
    sres@tmYLMinorOn = True
    sres@tmYLMinorValues = ispan(-40,40,4)*0.01
    sres@tmYLMode = "Explicit"
    plot(1) = gsn_csm_xy(wks,ispan(1, 6, 1),number(1,:) - number(0,:),sres)   ; draw each time series

    ; plot(1) = gsn_csm_xy(wks,ispan(1, 6, 1),number(1,:),sres)   ; draw each time series

    ; lbres                    = True          ; labelbar only resources
    ; lbres@vpWidthF           = 0.3           ; labelbar width
    ; lbres@vpHeightF          = 0.1           ; labelbar height
    ; lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
    ; lbres@lbFillColors       = (/"red","black"/)
    ; lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    ; lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
    ; lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    ; lbres@lbPerimOn          = False
    ; lbres@lgPerimColor    = "white"
    ; labels = (/"MSE","CpT"/)
    ; gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.23,lbres) ; draw right labelbar column

    ; lbres@lbFillColors       := (/"green","blue"/)
    ; labels := (/"Lvq","gz"/)
    ; gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.23,lbres) ; draw left labelbar column

    frame(wks)      

end


