;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
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

begin
season_r = (/"spring","summer","autumn","winter"/)
season = (/"MAM","JJA","SON","DJF"/)
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h1.0101-4012.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarimg_" + season_r +"_notest/"+"F_2000_atl_"+season_r+"_notest.cam.h1.0101-4012.nc" 
  
  ; wks=gsn_open_wks("pdf", "CESM-F2000-u,v,OMEGA-alt-warming-foul-season-global")
  
  ; plot = new(4,graphic)
  ; iplot = new(4,graphic)
  ; str3 = ispan(1,4,1) + ""


do i = 0,3
  f1=addfile(fils1, "r")
  f2=addfile(fils2(i), "r")
  print("读入file"+fils2(i))

  date=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.601.and.date.le.4012)

  ;;;;U,V
  u_CT  =  f1->U(timeind,{992},:,:)
  v_CT  =  f1->V(timeind,{992},:,:)
  
  u_id  =  f2->U(timeind,{992},:,:)
  v_id  =  f2->V(timeind,{992},:,:)

  
  dim  = dimsizes(u_id)
  
  u_diff = u_CT(0:11,:,:)


  v_diff = u_diff

  printVarSummary(u_diff)

 

  m = (/2,5,8,0/)

  if season(i).ne."DJF" then
    do j = 0,2
      u_diff(j+m(i),:,:) = dim_avg_n_Wrap(u_id(j+m(i):dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(j+m(i):dim(0)-1:12,:,:),0)
      v_diff(j+m(i),:,:) = dim_avg_n_Wrap(v_id(j+m(i):dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(j+m(i):dim(0)-1:12,:,:),0)
    end do 
  else 
      u_diff(11,:,:) = dim_avg_n_Wrap(u_id(11:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(11:dim(0)-1:12,:,:),0)
      u_diff(0,:,:) = dim_avg_n_Wrap(u_id(0:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(0:dim(0)-1:12,:,:),0)
      u_diff(1,:,:) = dim_avg_n_Wrap(u_id(1:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(u_CT(1:dim(0)-1:12,:,:),0)

      v_diff(11,:,:) = dim_avg_n_Wrap(v_id(11:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(11:dim(0)-1:12,:,:),0)
      v_diff(0,:,:) = dim_avg_n_Wrap(v_id(0:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(0:dim(0)-1:12,:,:),0)
      v_diff(1,:,:) = dim_avg_n_Wrap(v_id(1:dim(0)-1:12,:,:),0) - dim_avg_n_Wrap(v_CT(1:dim(0)-1:12,:,:),0)
  end if 


    ; u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season(i)), 0)
    ; u_id_clm = dim_avg_n_Wrap(month_to_season(u_id, season(i)), 0)

    ; u_diff = u_CT_clm
    ; u_diff = (/u_id_clm - u_CT_clm/)

    ; v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, season(i)), 0)
    ; v_id_clm = dim_avg_n_Wrap(month_to_season(v_id, season(i)), 0)

    ; v_diff = v_CT_clm
    ; v_diff = (/v_id_clm - v_CT_clm/)
    u_diff@_FillValue = 10^9
    v_diff@_FillValue = 10^9
;;;;;

end do 

;;-----------------------写入nc文件---------------------------
    
    ; system("/bin/rm -f /home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_"+season(i)+".nc")   ; remove any pre-existing file
    ; ncdf = addfile("/home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind_"+season(i)+".nc","c")  ; open output netCDF file


    system("/bin/rm -f /home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/Project/SCS-rain/POP_model_experiment/cam_atlwarming_wind.nc","c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->u_diff  =  u_diff                           ; 4D               
       ncdf->v_diff  = v_diff                       ; 3D
       
end 
