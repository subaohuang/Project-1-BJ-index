function area_avg(data,long_name)
local data_a
begin 
if (data@long_name.eq."Sea surface temperature") then
data_a=data@scale_factor*data+data@add_offset-274.15
else 
data_a=data@scale_factor*data+data@add_offset
end if
data_a=doubletofloat(-data_a)
copy_VarMeta(data,data_a )
data_a@long_name=long_name
data_a_new = data_a(latitude|:,longitude|:,time|:)
data_a_new&time=ispan(1,432,1)

data_avg_area=dim_avg_n_Wrap(data_a_new, (/0,1/))
data_avg_area2=doubletofloat(data_avg_area)
copy_VarMeta(data_a_new(1,1,:), data_avg_area2)

return(data_avg_area2)
end 

;------------------------------------------------------------------
;------------------------------------------------------------------
begin 
;读入数据

diri="/media/sf_share/era/data/"
fils1=systemfunc("ls "+diri+"era_u_sst_slp.nc")
f1=addfile(fils1, "r")
fils2=systemfunc("ls "+diri+"era_w_rq.nc")
f2=addfile(fils2, "r")
sst=f1->sst(0:431,40:80,:)
msl=f1->msl(0:431,40:80,:)
u10=f1->u10(0:431,40:80,:)
rq=f2->r(0:431,0,40:80,:)

wks2=gsn_open_wks("png", "/media/sf_share/era/era_u10_series")

sst_avg_area=area_avg(u10,"u10_avg_area")

sst_rc=regCoef(sst_avg_area&time,sst_avg_area)

data=new((/3,dimsizes(sst_avg_area&time)/),float)

data(0,:) = sst_avg_area
data(1,:) = sst_rc*sst_avg_area&time+sst_rc@yintercept
data(2,:) =avg(sst_avg_area)

;---------------------------------------------------
;画图 
;-------------------------------------------------------

resxy=True
resxy@gsnDraw=False
resxy@gsnFrame=False
resxy@tmXBMode="Explicit"
resxy@trXMaxF=456
resxy@tmXBValues=ispan(13,432,60)
resxy@tmXBLabels=ispan(1980,2015,5)

resxy@xyLineThicknesses = (/ 3.0,2.0,2.0/)          ; make second line thicker
resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
resxy@xyDashPatterns    =(/0,0,2/)


resxy@vpHeightF=0.3
resxy@vpWidthF=0.7
resxy@vpYF=0.9
resxy@vpXF=0.15

resxy@tmXTOn=False
resxy@tmYROn=False
resxy@tmYLMinorOn=False


resxy@tiYAxisString = "-u/m*s~S~-1~N"



; resxy@tmLabelAutoStride = True                    ; nice stride on labels
plot2=gsn_csm_xy(wks2, sst_avg_area&time,data,resxy)


draw(plot2)
frame(wks2)
end
