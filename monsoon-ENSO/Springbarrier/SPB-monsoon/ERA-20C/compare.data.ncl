  begin 
year = 1900
    fils1 := "/home/yangsong3/data-observation/ERA-20C/pressure/uwnd/uwnd.20C."+year+".nc"
    fils2 := "/home/yangsong3/data-observation/ERA-20C/pressure/uwnd/test.nc"
    
    f1 = addfile(fils1, "r")
    f2 = addfile(fils2, "r")
   
;;;;U,V
    u_200  = short2flt(f1->u(:,{850},:,:))
    u_200_double  = f2->u(:,{850},:,:)

print(u_200(4,30,30))
print(u_200_double(4,30,30))
end 
;;;
