; begin
;    grib_in= addfile("/mnt/j/JRA-25/anl_p25.033_ugrd.195801_195812.grb","r")
;    names  = getfilevarnames(grib_in)   ; extract ALL variable names
;    print(names)
;    wund=grib_in->UGRD_GDS0_ISBL_S123
;    printVarSummary(wund)
; end
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
;
; Create array to contain dataitem ids.
;
  ncurve  = 10
  npts    = 100
  dataids = new(ncurve-1,graphic)
  addids  = new(ncurve-1,graphic)

;
; Initialize some dummy data for the xy plot.
; 

  ii = tofloat(conform_dims((/ncurve,npts/),ispan(0,npts-1,1),1))
  jj = tofloat(conform_dims((/ncurve,npts/),ispan(1,ncurve,1),0))

  pi = 3.14159
  y = jj*sin((2.*ii*pi)/((npts-((jj-1)*10)-1)))
  delete([/ii,jj/])

;
; Begin graphics section.
;
  wks = gsn_open_wks("x11","adddata")    ; open workstation

  res                  = True
  res@gsnMaximize      = True
  res@gsnDraw          = False
  res@gsnFrame         = False

  res@trYMinF          = min(y)    ; Make sure we have room to plot
  res@trYMaxF          = max(y)    ; all curves.

  res@xyLineColor      = 2
  res@xyLineThicknessF = 2.

;
; Create XY plot, but don't draw it yet.
;
  xy = gsn_y(wks,y(0,:),res)
;
; Loop through each curve. For each remaining curve, create
; the data object and add it to the existing XY plot.
;
  do i=1,ncurve-1
;
; Create data object.
;
    dataids(i-1) = create "xyData"+i coordArraysClass defaultapp
      "caYArray": y(i,0:(npts-i*10)-1)
    end create
;
; Add to existing plot.
;
    addids(i-1) = NhlAddData(xy,"xyCoordData",dataids(i-1))
;
; Set some resources for this particular curve.
;
    setvalues addids(i-1)
      "xyLineColor"      : i+3
      "xyLineThicknessF" : 1.02*(i+2)
    end setvalues
  end do

;
; Now draw the plot and advance the frame. You should see all
; ten curves, each with a different color and thickness.
; 
  draw(xy)
  frame(wks)

;
; Remove every other data item and redraw.
; You should see only five curves now.
;
  NhlRemoveData(xy,"xyCoordData",dataids(::2))
  draw(xy)
  frame(wks)

end