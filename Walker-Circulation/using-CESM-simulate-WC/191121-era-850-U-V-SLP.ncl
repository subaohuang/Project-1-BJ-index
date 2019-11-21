begin
; fils1  =  "/home/ys17-19/lsh/data/era/era_pre_global.nc"
fils1  =  "/home/ys17-19/lsh/data/era/era_sur_global.nc"


    f1 = addfile(fils1, "r")
    

    date  =  cd_calendar(f1->time, 1)

    timeind  =  ind(date.le.201212.and.date.ge.197901)
   

    u_heat = short2flt(f1->u10(timeind,:,:))

    v_heat = short2flt(f1->v10(timeind,:,:))
   

    season = (/"MAM","JJA","SON","DJF"/)
    
    dim = dimsizes(u_heat)

    u_heat_mon = new((/12,34,dim(1),dim(2)/),typeof(u_heat))
    
    do i = 0,11,1
    u_heat_mon(i,:,:,:) = u_heat_mon(i:dim(0)-1:12,:,:)
    end do
printVarSummary(u_heat_mon)
exit()

    do i = 0,3   
    u_heat_mam = month_to_season(u_heat, season(i))

    v_heat_mam = month_to_season(v_heat, season(i))


    u_heat_mam_pa = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{150:210}), (/1,2/))

    u_heat_mam_id = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{70:100}), (/1,2/))

    corr_heat = escorc_n(u_heat_mam_pa, u_heat_mam_id, 0, 0)


    print("corr_heat="+corr_heat)
    end do 


end    





;;;
