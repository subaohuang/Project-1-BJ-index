
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"
begin
  

; =================================================================
; Access the area of interest via coordinate subscripting
; Note that by using eofunc_n_Wrap and eofunc_ts_n_Wrap, we 
; don't have to reorder data first.
; =================================================================



  ; ======================================
  filename = "xd.csv";
  values_1d = asciiread(filename,-1,"float") ;将数据按照一维变量的格式读取，变量类型为字符串
  ncols = dimsizes(str_split(values_1d(0),",")) ;先将数据每一行用逗号分隔，再读取它的维数 
  nrows = dimsizes(values_1d) ;将读取数据维数 ;---根据需要把数据重新排列为2维数据，并且转换为整
  values_2d = onedtond(values_1d,(/nrows,ncols/)) ;先把一维变量转换为多维变量，再强制转换为整型数据
  ;write_matrix(values_2d,ncols+"I5",False) ;
printVarSummary(values_2d)
exit()
  values_2d@_FillValue = 999;
  ;printVarSummary(values_2d);
  ;values_2d(1,1) = 100;
  print(values_2d(1,1));

  ;time = values_2d(:,1);








  time = (/15769752,17645640/);

  time@units = "hours since 1-1-1 00:00:0.0" ; 设定时间单位
  utc_date = cd_calendar(time, 0)  ; 转换到世界时
  year   = tointeger(utc_date(:,0))
  month  = tointeger(utc_date(:,1))
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  minute = tointeger(utc_date(:,4))
  second = utc_date(:,5)
  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)
  date_str = sprinti("%0.2i ", hour) + sprinti("%0.2i ", day) + \
             month_abbr(month) + " "  + sprinti("%0.4i", year)
  ;print(date_str)









  ;xw     = wSLP(time|:,{lat|latS:latN},{lon|lonL:lonR})

  ;eof    = eofunc_n_Wrap(xw, 3, True, 0)
  ;eof_ts = eofunc_ts_n_Wrap (xw, eof, False, 0)

  ;printVarSummary( eof )                         ; examine EOF variables
  ;printVarSummary( eof_ts )
end