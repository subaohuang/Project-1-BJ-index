;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd
;;linshuheng 2019.12.24
;**********************************************************
begin
  ;;;read era data
		path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
		var_name = (/"uwnd","vwnd","shum"/)
		filein = path_wind + var_name +".mon.mean.nc"

		f1 = addfile(filein(0),"r")

		date=cd_calendar(f1->time, 1)
		timeind=ind(date.ge.197901.and.date.le.201312)

		u = short2flt(f1->u(timeind,:,:,:))

		f2 = addfile(filein(1),"r")
		v = short2flt(f2->v(timeind,:,:,:))

		
		f3 = addfile(filein(2),"r")
		date :=cd_calendar(f3->time, 1)
		timeind :=ind(date.ge.197901.and.date.le.201312)
		q = short2flt(f3->q(timeind,:,:,:))
          

   ;;;


   ;;;

   






