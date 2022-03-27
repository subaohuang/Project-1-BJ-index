;;*************************************************
;; main program 
;;*************************************************************
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

begin
    startmon = 101
    endmon = 99712
    startyear = 1
    endyear = 997

    lat1 = -30
    lat2 = 30
    lon1 = -180
    lon2 = 180

; ; ;;;读入海温数据
    f1 = addfile("ssta_picontrol.nc", "r")
    ssta = f1->ssta
;;;;;;

; ; ;;;;
    f1 = addfile("ssta_picontrol_highpass_bw.nc", "r")
    ssta_highpass = f1->ssta
;;;;

printVarSummary(ssta)
printVarSummary(ssta_highpass)
   nino3_4 = dim_avg_n_Wrap(ssta(:,{-5:5},{190:240}), (/1,2/))
   nino3_4_highpass = dim_avg_n_Wrap(ssta_highpass(:,{-5:5},{190:240}), (/1,2/))

   nino3_4_jjas = month_to_season(nino3_4, "JJA")
   nino3_4_highpass_jjas = month_to_season(nino3_4_highpass, "JJA")

x = ispan(1,dimsizes(nino3_4_jjas),1)
wks = gsn_open_wks("pdf", "nino3_4_jjas_index")
res = True
res@xyLineColors = (/"black","red"/)
plot = gsn_csm_xy(wks, x, (/nino3_4_jjas,nino3_4_highpass_jjas/), res)

end 