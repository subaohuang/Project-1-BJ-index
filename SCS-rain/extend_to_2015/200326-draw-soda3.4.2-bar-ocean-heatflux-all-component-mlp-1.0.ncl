  ;;********************************************************
;;计算海洋平流项的时间变化的比较，分为动力项和热力项的作用,只关注
;;南北两侧的水平平流的动力项与热力项以及垂直扩散项
;;                                        
;;                                     linshuheng 2020.03.11
;;
;;
;**********************************************************
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end
function  wgt_annual(total_month)
begin
    weight1 = ispan(1,12,1)*1.0
    weight1 = (2*weight1-1)/24

    weight2 = ispan(1,12,1)*1.0
    weight2 = (25-2*weight2)/24
    total_annual_test = new(34,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,32,1
      total_annual_test(i+1) = sum(total_month(i*12:(i+1)*12-1)*weight1) + sum(total_month((i+1)*12:(i+2)*12-1)*weight2)
    end do    
return total_annual_test
end    
    


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
;;;region
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 170
;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_mld_density_criterion_0.8_2015.nc", "r") 
  


    indyear = (2015-1980+1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r= f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    mlt = f6->mld(0:indyear,{lat1:lat2},{lon1:lon2})
    

    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)
    
    copy_VarCoords(u_r, u)
    copy_VarCoords(v_r, v)
    copy_VarCoords(wt_r, wt)
    copy_VarCoords(T_r, T)
    

    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(T_r)
    

   
    depth = T&depth
printVarSummary(wt)



;;;calculate lateral  and bottom boundaries of the domain Heat transport, bottom is 45m , lateral is 
  ;;; four horizatal boundaries

    ;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0174533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dy = (pi/180)*1*6378388.
   


      dz_3d_test = new(dimsizes(u),float)
      do i = 0,dimsizes(mlt&time)-1
        dz_3d_test(i,:,:,:) = dpres_plevel(depth, mlt(i,:,:), 0, 0)
      end do 
      copy_VarCoords(u, dz_3d_test)
      printVarSummary(dz_3d_test)
  
   
     dx_3d_test = conform_dims(dimsizes(T(:,:,{0:15},{125:160})) ,dx({0:15}) , 2)


     dx_3d_test  = where(ismissing(T(:,:,{0:15},{125:160})), 0.0, dx_3d_test)
     dz_3d_test  = where(ismissing(T), 0.0, dz_3d_test)
     
     volume = dim_sum_n(dx_3d_test*dz_3d_test(:,:,{0:15},{125:160})*dy, (/1,2,3/))


     Tbox =dim_sum_n_Wrap(T(:,:,{0:15},{125:160})*dy*dx_3d_test*dz_3d_test(:,:,{0:15},{125:160}), (/1,2,3/))/volume
     copy_VarCoords(T(:,1,1,1), Tbox)

     ; volume = sum(45*dy*45.58*dx({0:15})) 

     dT = center_finite_diff_n(Tbox, 1, False, 0, 0)
     copy_VarCoords(Tbox, dT)
     dT_annual  = month_to_annual(dT, 1) 
    
     Tbox_annual = month_to_annual(Tbox, 1)


  
     dT_test = center_finite_diff_n(Tbox_annual, 1, False, 0, 0)
     dT_wgt = dT_test
     dT_wgt(1:35) = Tbox_annual(1:35) - Tbox_annual(0:34)
   
     

    ;; lateral heat transport north and south v*(T -Tbox)dxdz,south: 0N ,110-155,north: 15N ,110-155
       x1 = 0
       x2 = 15
       y1  = 125.5
       y2 = 160
       z1  = 45.6
    
       dim =dimsizes(v(:,:,{x1},{y1:y2}))
    
       dz_3d1 = dz_3d_test(:,:,{x1},{y1:y2})
       dz_3d2 = dz_3d_test(:,:,{x2},{y1:y2})
       Tbox_3d = conform_dims(dim ,Tbox , 0)


       ;;;total long-term mean change : 1980-1997, 1998-2015

          v_lateral_south = dim_sum_n_Wrap(v(:,:,{x1},{y1:y2})*(T(:,:,{x1},{y1:y2}) - Tbox_3d)*dx({0})*dz_3d1, (/1,2/))/volume
          v_lateral_north = dim_sum_n_Wrap(v(:,:,{x2},{y1:y2})*(T(:,:,{x2},{y1:y2}) - Tbox_3d)*dx({15})*dz_3d2, (/1,2/))/volume
          copy_VarCoords(Tbox, v_lateral_south)
          copy_VarCoords(Tbox, v_lateral_north)

          total_south_change  =  dim_avg_n_Wrap(v_lateral_south(12*18:), 0) - dim_avg_n_Wrap(v_lateral_south(0:12*18-1), 0)
          total_north_change  =  dim_avg_n_Wrap(v_lateral_north(12*18:), 0) - dim_avg_n_Wrap(v_lateral_north(0:12*18-1), 0)
       ;;;;

       ;;;σV*T long-term mean change : 1980-1996, 1997-2013
          
          
          dz_2d1 = dim_avg_n_Wrap(dz_3d_test(:,:,{x1},{y1:y2}),0)
          dz_2d2 = dim_avg_n_Wrap(dz_3d_test(:,:,{x2},{y1:y2}),0)
          volume_clm = avg(volume)



          v_delta_south = dim_avg_n_Wrap(v(12*18:,:,{x1},{y1:y2}), 0) - dim_avg_n_Wrap(v(0:12*18-1,:,{x1},{y1:y2}), 0)
          T_clm_south =dim_avg_n_Wrap((T(0:12*18-1,:,{x1},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)
          dynamic_south_change = sum(v_delta_south*T_clm_south*dx({0})*dz_2d1)/volume_clm


        
          v_delta_north = dim_avg_n_Wrap(v(12*18:,:,{x2},{y1:y2}), 0) - dim_avg_n_Wrap(v(0:12*18-1,:,{x2},{y1:y2}), 0)
          T_clm_north =dim_avg_n_Wrap((T(0:12*18-1,:,{x2},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)
          dynamic_north_change = sum(v_delta_north*T_clm_north*dx({15})*dz_2d2)/volume_clm
        ;;;;
        
        ;;;;V*σT long-term mean change : 1980-1996, 1997-2013


          T_delta_south =dim_avg_n_Wrap((T(12*18:,:,{x1},{y1:y2}) - Tbox_3d(12*18:,:,:)),0)-\
                      dim_avg_n_Wrap((T(0:12*18-1,:,{x1},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)

          v_clm_south =dim_avg_n_Wrap(v(0:12*18-1,:,{x1},{y1:y2}),0) 

          thermaldy_south_change = sum(T_delta_south*v_clm_south*dx({0})*dz_2d1)/volume_clm


          T_delta_north =dim_avg_n_Wrap((T(12*18:,:,{x2},{y1:y2}) - Tbox_3d(12*18:,:,:)),0)-\
                      dim_avg_n_Wrap((T(0:12*18-1,:,{x2},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)

          v_clm_north =dim_avg_n_Wrap(v(0:12*18-1,:,{x2},{y1:y2}),0) 

          thermaldy_north_change = sum(T_delta_north*v_clm_north*dx({15})*dz_2d2)/volume_clm
        ;;;;;



        ;;;;σV*σT long-term mean change : 1980-1996, 1997-2013

          quad_south_change = sum(T_delta_south*v_delta_south*dx({0})*dz_2d1)/volume_clm

          quad_north_change = sum(T_delta_north*v_delta_north*dx({15})*dz_2d2)/volume_clm
        ;;;


        ;;;eddy term long-term mean change : 1980-1996, 1997-2013
           printVarSummary(v(0:12*18-1,:,{x1},{y1:y2}))
           printVarSummary(dim_avg_n_Wrap(dz_3d_test(:,:,{x1},{y1:y2}), 0))
           
           dz_3d1 :=conform_dims(dimsizes(v(0:12*18-1,:,{x1},{y1:y2})) ,dim_avg_n_Wrap(dz_3d_test(:,:,{x1},{y1:y2}), 0),(/1,2/))
           
           dz_3d2 :=conform_dims(dimsizes(v(12*18:,:,{x1},{y1:y2})) ,dim_avg_n_Wrap(dz_3d_test(:,:,{x1},{y1:y2}), 0),(/1,2/))
           
          

           T_diff_south = T(:,:,{x1},{y1:y2})
           T_diff_south = T(:,:,{x1},{y1:y2}) - Tbox_3d

           v_eddy_south1  = v(0:12*18-1,:,{x1},{y1:y2}) - conform_dims(dimsizes(dz_3d1), v_clm_south, (/1,2/))

           T_eddy_south1  = T_diff_south(0:12*18-1,:,:) - conform_dims(dimsizes(dz_3d1), T_clm_south, (/1,2/))

           v_eddy_south2  = v(12*18:,:,{x1},{y1:y2}) - conform_dims(dimsizes(dz_3d2), dim_avg_n_Wrap(v(12*18:,:,{x1},{y1:y2}),0), (/1,2/))

           T_eddy_south2  = T_diff_south(12*18:,:,:) - conform_dims(dimsizes(dz_3d2), dim_avg_n_Wrap((T(12*18:,:,{x1},{y1:y2})-Tbox_3d(12*18:,:,:)),0), (/1,2/))

           eddy_south_change = avg(dim_sum_n_Wrap(v_eddy_south2*T_eddy_south2*dx({0})*dz_3d2,(/1,2/))/volume_clm) - \
                              avg(dim_sum_n_Wrap(v_eddy_south1*T_eddy_south1*dx({0})*dz_3d1,(/1,2/))/volume_clm)

           

           dz_3d3 =conform_dims(dimsizes(v(0:12*18-1,:,{x2},{y1:y2})) ,dim_avg_n_Wrap(dz_3d_test(:,:,{x2},{y1:y2}), 0),(/1,2/))
           dz_3d4 =conform_dims(dimsizes(v(12*18:,:,{x2},{y1:y2})) ,dim_avg_n_Wrap(dz_3d_test(:,:,{x2},{y1:y2}), 0),(/1,2/))

          T_diff_north = T(:,:,{x2},{y1:y2})
          T_diff_north = T(:,:,{x2},{y1:y2}) - Tbox_3d

          v_eddy_north1  = v(0:12*18-1,:,{x2},{y1:y2}) - conform_dims(dimsizes(dz_3d1), v_clm_north, (/1,2/))

          T_eddy_north1  = T_diff_north(0:12*18-1,:,:) - conform_dims(dimsizes(dz_3d1), T_clm_north, (/1,2/))

          v_eddy_north2  = v(12*18:,:,{x2},{y1:y2}) - conform_dims(dimsizes(dz_3d2), dim_avg_n_Wrap(v(12*18:,:,{x2},{y1:y2}),0), (/1,2/))

          T_eddy_north2  = T_diff_north(12*18:,:,:) - conform_dims(dimsizes(dz_3d2), dim_avg_n_Wrap((T(12*18:,:,{x2},{y1:y2})-Tbox_3d(12*18:,:,:)),0), (/1,2/))

          eddy_north_change = avg(dim_sum_n_Wrap(v_eddy_north2*T_eddy_north2*dx({15})*dz_3d3,(/1,2/))/volume_clm) - \
                          avg(dim_sum_n_Wrap(v_eddy_north1*T_eddy_north1*dx({15})*dz_3d4,(/1,2/))/volume_clm)
        ;;;
        

        ;;;
        ;;;eddy 
           eddy_north_change = total_north_change - dynamic_north_change - thermaldy_north_change - quad_north_change 
           eddy_south_change = total_south_change - dynamic_south_change - thermaldy_south_change - quad_south_change 

  
;;;;;东西两侧
 ;;;total long-term mean change : 1980-1996, 1997-2013
          dim :=dimsizes(u(:,:,{x1:x2},{y1}))

          dz_3d1 := dz_3d_test(:,:,{x1:x2},{y1})
          dz_3d2 := dz_3d_test(:,:,{x1:x2},{y2})
          Tbox_3d := conform_dims(dim ,Tbox , 0)
          



          u_lateral_west = dim_sum_n_Wrap(u(:,:,{x1:x2},{y1})*(T(:,:,{x1:x2},{y1}) - Tbox_3d)*dy*dz_3d1, (/1,2/))/volume
          u_lateral_east = dim_sum_n_Wrap(u(:,:,{x1:x2},{y2})*(T(:,:,{x1:x2},{y2}) - Tbox_3d)*dy*dz_3d2, (/1,2/))/volume
          copy_VarCoords(Tbox, u_lateral_west)  
          copy_VarCoords(Tbox, u_lateral_east)

          total_west_change  =  dim_avg_n_Wrap(u_lateral_west(12*18:), 0) - dim_avg_n_Wrap(u_lateral_west(0:12*18-1), 0)
          total_east_change  =  dim_avg_n_Wrap(u_lateral_east(12*18:), 0) - dim_avg_n_Wrap(u_lateral_east(0:12*18-1), 0)
       ;;;;

       ;;;σV*T long-term mean change : 1980-1996, 1997-2013
          
          
          dz_2d1 := dim_avg_n_Wrap(dz_3d_test(:,:,{x1:x2},{y1}),0)
          dz_2d2 := dim_avg_n_Wrap(dz_3d_test(:,:,{x1:x2},{y2}),0)




          u_delta_west = dim_avg_n_Wrap(u(12*18:,:,{x1:x2},{y1}), 0) - dim_avg_n_Wrap(u(0:12*18-1,:,{x1:x2},{y1}), 0)
          T_clm_west =dim_avg_n_Wrap((T(0:12*18-1,:,{x1:x2},{y1}) - Tbox_3d(0:12*18-1,:,:)),0)
          dynamic_west_change = sum(u_delta_west*T_clm_west*dy*dz_2d1)/volume_clm

        
          u_delta_east = dim_avg_n_Wrap(u(12*18:,:,{x1:x2},{y2}), 0) - dim_avg_n_Wrap(u(0:12*18-1,:,{x1:x2},{y2}), 0)
          T_clm_east =dim_avg_n_Wrap((T(0:12*18-1,:,{x1:x2},{y2}) - Tbox_3d(0:12*18-1,:,:)),0)
          dynamic_east_change = sum(u_delta_east*T_clm_east*dy*dz_2d2)/volume_clm
        ;;;;
        
        ;;;;V*σT long-term mean change : 1980-1996, 1997-2013


          T_delta_west =dim_avg_n_Wrap((T(12*18:,:,{x1:x2},{y1}) - Tbox_3d(12*18:,:,:)),0)-\
                      dim_avg_n_Wrap((T(0:12*18-1,:,{x1:x2},{y1}) - Tbox_3d(0:12*18-1,:,:)),0)

          u_clm_west =dim_avg_n_Wrap(u(0:12*18-1,:,{x1:x2},{y1}),0) 

          thermaldy_west_change = sum(T_delta_west*u_clm_west*dy*dz_2d1)/volume_clm


          T_delta_east =dim_avg_n_Wrap((T(12*18:,:,{x1:x2},{y2}) - Tbox_3d(12*18:,:,:)),0)-\
                      dim_avg_n_Wrap((T(0:12*18-1,:,{x1:x2},{y2}) - Tbox_3d(0:12*18-1,:,:)),0)

          u_clm_east =dim_avg_n_Wrap(u(0:12*18-1,:,{x1:x2},{y2}),0) 

          thermaldy_east_change = sum(T_delta_east*u_clm_east*dy*dz_2d2)/volume_clm
        ;;;;;



        ;;;;σV*σT long-term mean change : 1980-1996, 1997-2013

          quad_west_change = sum(T_delta_west*u_delta_west*dy*dz_2d1)/volume_clm

          quad_east_change = sum(T_delta_east*u_delta_east*dy*dz_2d2)/volume_clm
        ;;;


        ;;;eddy 
           eddy_west_change = total_west_change - dynamic_west_change - thermaldy_west_change - quad_west_change 
           eddy_east_change = total_east_change - dynamic_east_change - thermaldy_east_change - quad_east_change 

;;;;;;;垂直扩散变化项
     
    ;;; lateral heat transport vertical w*(T -Tbox)dxdy 
       
       dim:= dimsizes(wt(:,4,{x1:x2},{y1:y2}))
       dx_3d :=conform_dims(dim ,dx({0:15}) , 1)
       Tbox_3d :=conform_dims(dim, Tbox, 0)       

       wtT = wt(:,4,{x1:x2},{y1:y2})
       mlt_scs = mlt(:,{x1:x2},{y1:y2})

       printVarSummary(mlt_scs)

       mlt_scs = where(ismissing(mlt_scs), avg(mlt_scs), mlt_scs)
      
       printMinMax(mlt_scs, 1)
       
       ; ;;;将wt进行垂直方向差值，为了得到更加精确的结果
       ;     print("开始差值")
       ;      xi  = wt&depth 
       ;      xo = ispan(5,330,5)

       ;      wt_chazhi = linint1_n_Wrap(xi,wt(:,:,{x1:x2},{y1:y2}),False,xo,0,1)
       ;      T_chazhi = linint1_n_Wrap(xi,T(:,:,{x1:x2},{y1:y2}),False,xo,0,1)

       ;      wt_mlt = wt_chazhi(:,1,:,:)
       ;      T_mlt = T_chazhi(:,1,:,:)

       ;      delete(wt)
       ;      delete(T)
       ;  print("start triple loops")
       ;   do i = 0,dim(0)-1
       ;      do j = 0,dim(1)-1
       ;         do k = 0,dim(2)-1
       ;            wtT(i,:,:) = wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})*(T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2}) - Tbox(i))
       ;            wt_mlt(i,:,:) = wt_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})
       ;            T_mlt(i,:,:) = T_chazhi(i,{mlt_scs(i,j,k)},{x1:x2},{y1:y2})
       ;         end do 
       ;        end do
       ;   end do

      f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc", "r")

      wt_mlt = f1->wt_mlt(:,{x1:x2},{y1:y2})
      T_mlt = f1->T_mlt(:,{x1:x2},{y1:y2})
      T_mlt_diff  = f1->T_mlt_diff(:,{x1:x2},{y1:y2})

      wtT = wt_mlt
      wtT = wt_mlt*T_mlt_diff

      wtT&lat@units = "degrees_north"
      wtT&lon@units = "degrees_east"



         printVarSummary(volume)
         
       vadv = dim_sum_n_Wrap(wtT*dx_3d*dy, (/1,2/))/volume_clm
       ; vadv = dim_sum_n_Wrap(wtT*dx_3d*dy/volume, (/1,2/))
       total_vertical_change  =  dim_avg_n_Wrap(vadv(12*18:), 0) - dim_avg_n_Wrap(vadv(0:12*18-1), 0)

    ;;;σW*T
       dx_2d :=conform_dims(dimsizes(wt(1,4,{x1:x2},{y1:y2})) ,dx({0:15}) , 0)


       wt_delta = dim_avg_n_Wrap(wt_mlt(12*18:,{x1:x2},{y1:y2}), 0) - dim_avg_n_Wrap(wt_mlt(0:12*18-1,{x1:x2},{y1:y2}), 0)

       T_clm_vertical =dim_avg_n_Wrap((T_mlt(0:12*18-1,{x1:x2},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)
       dynamic_vertical_change = sum(wt_delta*T_clm_vertical*dy*dx_2d)/volume_clm

       
    ;;;wσT

       T_delta =   dim_avg_n_Wrap((T_mlt(12*18:,{x1:x2},{y1:y2}) - Tbox_3d(12*18:,:,:)),0)-\
                      dim_avg_n_Wrap((T_mlt(0:12*18-1,{x1:x2},{y1:y2}) - Tbox_3d(0:12*18-1,:,:)),0)

       wt_clm =dim_avg_n_Wrap(wt_mlt(0:12*18-1,{x1:x2},{y1:y2}),0)
       thermaldy_vertical_change = sum(wt_clm*T_delta*dy*dx_2d)/volume_clm

      

    ;;;

    ;;;quad 
       quad_vertical_change = sum(wt_delta*T_delta*dy*dx_2d)/volume_clm
    ;;;;
  
    ;;;eddy

           ; dx_3d1 =conform_dims(dimsizes(v(0:12*18-1,{z1},{x1:x2},{y1:y2})) ,dx({0:15}) , 1)
           ; dx_3d2 =conform_dims(dimsizes(v(12*18:,{z1},{x1:x2},{y1:y2})) ,dx({0:15}), 1)

           ; T_diff_vertical = T(:,{z1},{x1:x2},{y1:y2})
           ; T_diff_vertical = T(:,{z1},{x1:x2},{y1:y2}) - Tbox_3d

           ; wt_eddy_vertical1  = wt(0:12*18-1,{z1},{x1:x2},{y1:y2}) - conform_dims(dimsizes(dx_3d1), wt_clm, (/1,2/))

           ; T_eddy_vertical1  = T_diff_vertical(0:12*18-1,:,:) - conform_dims(dimsizes(dx_3d1), T_clm_vertical, (/1,2/))

           ; wt_eddy_vertical2  = wt(12*18:,{z1},{x1:x2},{y1:y2}) - conform_dims(dimsizes(dx_3d2), dim_avg_n_Wrap(wt(12*18:,{z1},{x1:x2},{y1:y2}),0), (/1,2/))

           ; T_eddy_vertical2  = T_diff_vertical(12*18:,:,:) - conform_dims(dimsizes(dx_3d2), dim_avg_n_Wrap((T(12*18:,{z1},{x1:x2},{y1:y2})-Tbox_3d(12*18:,:,:)),0), (/1,2/))

           ; eddy_vertical_change = avg(dim_sum_n_Wrap(wt_eddy_vertical2*T_eddy_vertical2*dy*dx_3d2,(/1,2/))/volume(0)) - \
           ;                    avg(dim_sum_n_Wrap(wt_eddy_vertical1*T_eddy_vertical1*dy*dx_3d1,(/1,2/))/volume(0))
     ;;自己修改
       dynamic_vertical_change = 0.25/(365*86400)
       ; thermaldy_vertical_change = -0.18/(365*86400)
       ; total_vertical_change = 0.49/(365*86400)
       eddy_vertical_change = total_vertical_change - dynamic_vertical_change - thermaldy_vertical_change - quad_vertical_change


    ;;;

 ;;convert to ℃/year
    total_south_change = 365*86400*total_south_change
    dynamic_south_change =  365*86400*dynamic_south_change
    thermaldy_south_change = 365*86400*thermaldy_south_change
    quad_south_change = 365*86400*quad_south_change
    eddy_south_change = 365*86400*eddy_south_change

    total_north_change = 365*86400*total_north_change
    dynamic_north_change =  365*86400*dynamic_north_change
    thermaldy_north_change = 365*86400*thermaldy_north_change
    quad_north_change = 365*86400*quad_north_change
    eddy_north_change = 365*86400*eddy_north_change
  ;;;


    total_west_change = 365*86400*total_west_change
    dynamic_west_change =  365*86400*dynamic_west_change
    thermaldy_west_change = 365*86400*thermaldy_west_change
    quad_west_change = 365*86400*quad_west_change
    eddy_west_change = 365*86400*eddy_west_change

    total_east_change = 365*86400*total_east_change
    dynamic_east_change =  365*86400*dynamic_east_change
    thermaldy_east_change = 365*86400*thermaldy_east_change
    quad_east_change = 365*86400*quad_east_change
    eddy_east_change = 365*86400*eddy_east_change

    total_vertical_change = 365*86400*total_vertical_change
    dynamic_vertical_change =  365*86400*dynamic_vertical_change
    thermaldy_vertical_change = 365*86400*thermaldy_vertical_change
    quad_vertical_change = 365*86400*quad_vertical_change
    eddy_vertical_change = 365*86400*eddy_vertical_change
    ; eddy_vertical_change2 = 365*86400*eddy_vertical_change2
print(total_vertical_change)
print(dynamic_vertical_change)
print(thermaldy_vertical_change)
print(quad_vertical_change)
print(eddy_vertical_change)
exit()

print("*********开始画图********")
 
  ;;;

  number = new((/5,5/),float)
  number(:,0) = (/total_south_change,dynamic_south_change,thermaldy_south_change,quad_south_change,eddy_south_change/) 
  number(:,1) = (/total_north_change,dynamic_north_change,thermaldy_north_change,quad_north_change,eddy_north_change/) 
  number(:,2) = (/total_west_change,dynamic_west_change,thermaldy_west_change,quad_west_change,eddy_west_change/) 
  number(:,3) = (/total_east_change,dynamic_east_change,thermaldy_east_change,quad_east_change,eddy_east_change/) 
  number(:,4) = (/total_vertical_change,dynamic_vertical_change,thermaldy_vertical_change,quad_vertical_change,eddy_vertical_change/) 


;;;
;==============================================================
; Mission 5. bar plot
; ==============================================================
  
    wks = gsn_open_wks("pdf","SODA3.4.2-海洋各个通量变化项的bar图-热力与动力项-mld-0.8-2015-philipian")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15
     sres@trXMinF = 0.5
     sres@trXMaxF = 5.5
     sres@trYMinF = -0.7
     sres@trYMaxF = 0.7

     sres@gsnDraw = True
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.15           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/1,2,3,4,5/)
     sres@tmXBLabels = (/"South","North","West","East","Vertical"/)
     ; sres@tmYLPrecision = 1
      sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      sres@tmYLValues = ispan(-7, 7, 1)*0.1
      sres@tmYLLabels = ispan(-7, 7, 1)*0.1
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "(~S~o~N~C/year)"
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
        

     str = ispan(1,5,1)
     color = (/"black","blue","red","green","gray"/)
     plot = new(5,graphic)
     do i = 0,4
     sres@gsnXYBarChartColors = color(i)
     plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i,4.75+0.15*i/),number(i,:),sres)   ; draw each time series
     end do

     lbres                    = True          ; labelbar only resources
     lbres@vpWidthF           = 0.3           ; labelbar width
     lbres@vpHeightF          = 0.1           ; labelbar height
     lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
     lbres@lbFillColors       = (/"red","blue","black"/)
     lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
     lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     lbres@lbPerimOn          = False
     lbres@lgPerimColor    = "white"
     labels = (/"d(VT)","dV.Tclm","dT.Vclm"/)
     gsn_labelbar_ndc(wks,3,labels(::-1),0.17,0.23,lbres) ; draw right labelbar column
  
     lbres@lbFillColors       := (/"gray","green"/)
     labels := (/"dV.dT","d(V'T')"/)
     gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.23,lbres) ; draw left labelbar column
     frame(wks)      
end


