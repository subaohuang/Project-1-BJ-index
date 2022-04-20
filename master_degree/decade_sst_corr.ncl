begin
    diri="/mnt/f/share/era/data/"
    ipindex=asciiread("/mnt/f/share/era/data/ip_index_1950_2010.ascii", 732, "float")
   


    ipindex_annual=month_to_annual(ipindex, 1)
    ipindex_anu_std=dim_standardize_n(ipindex_annual,0,0)
    ipindex_anu_std!0="time"
    ipindex_anu_std&time=ispan(1,61,1)
     
    ipindex_anu_std_de=dtrend(ipindex_anu_std, True)
    copy_VarCoords(ipindex_anu_std, ipindex_anu_std_de)



    fils2=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc") 
    f2=addfile(fils2, "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.le.201012.and.date.ge.195001)
    sst=lonFlip(f2->sst(timeind,{20:-20},:))
    delete(sst&time)

    sst&time=date(timeind)
    
    sst_annual=month_to_annual(sst, 1)
    sst_annual&year=ispan(1, 61, 1)

    sst_de=dtrend_msg_n(ispan(0,dimsizes(sst_annual&year)-1,1), sst_annual, True, False, 0)
    copy_VarCoords(sst_annual, sst_de)

    sst_box1=dim_avg_n_Wrap(sst_annual(:,{10:-10},{90:150}), (/1,2/))   
    sst_box2=dim_avg_n_Wrap(sst_annual(:,{-4:-10},{240:270}), (/1,2/))
    ; sst_box3=dim_avg_n_Wrap(sst_annual(:,{2:-6},{60:90}), (/1,2/))
    sst_box3=dim_avg_n_Wrap(sst_annual(:,{2:-6},{40:80}), (/1,2/))
    sst_box4=dim_avg_n_Wrap(sst_annual(:,{5:-5},{140:170}), (/1,2/))
   
    sst_box5=dim_avg_n_Wrap(sst_annual(:,{5:-5},{140:170}), (/1,2/))   ;;整个印度洋洋盆，5到-5

    sst_box1_de=dim_avg_n_Wrap(sst_de(:,{10:-10},{80:170}), (/1,2/))
    sst_box2_de=dim_avg_n_Wrap(sst_de(:,{-4:-10},{240:270}), (/1,2/))
    ; sst_box3_de=dim_avg_n_Wrap(sst_de(:,{2:-6},{60:90}), (/1,2/))
    sst_box3_de=dim_avg_n_Wrap(sst_de(:,{2:-6},{40:80}), (/1,2/))
    sst_box4_de=dim_avg_n_Wrap(sst_de(:,{5:-5},{140:170}), (/1,2/))
    

    

    ipindex_anu_std_sm=runave_n_Wrap(ipindex_anu_std, 9, 0, 0)
    


    sst_box1_sm=runave_n_Wrap(sst_box1, 9, 0, 0)
    sst_box2_sm=runave_n_Wrap(sst_box2, 9, 0, 0)
    sst_box3_sm=runave_n_Wrap(sst_box3, 9, 0, 0)
    sst_box4_sm=runave_n_Wrap(sst_box4, 9, 0, 0)

    printVarSummary(sst_box1_sm)
    printVarSummary(ipindex_anu_std_sm)

    r1=escorc_n( ipindex_anu_std_sm, sst_box1_sm,0,0)
    
    r2=escorc_n( ipindex_anu_std_de, sst_box1_de,0,0)

    r3=escorc_n( ipindex_anu_std_sm, sst_box2_sm,0,0)

    r4=escorc_n( ipindex_anu_std_de, sst_box2_de,0,0)

    r5=escorc_n( ipindex_anu_std_sm, sst_box3_sm,0,0)

    r6=escorc_n( ipindex_anu_std_de, sst_box3_de,0,0)

    r7=escorc_n( ipindex_anu_std_de, sst_box4_de,0,0)

    r8=escorc_n( ipindex_anu_std_sm, sst_box4_sm,0,0)


    print("r1="+r1)
    
    print("r2="+r2)
     
    print("r3="+r3)

    print("r4="+r4)
   
    print("r5="+r5)

    print("r6="+r6)
    
    print("r7="+r7)
    
    print("r8="+r8)

;;;;超前滞后相关



end 


   

   













