;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
;************************************************
begin
  wks = gsn_open_wks("ncgm" ,"cyclones")             ; open a ps file
  wks = gsn_open_wks("pdf" ,"cyclones")             ; open a ps file
;  gsn_define_colormap(wks,"gsltod")             ; choose colormap
;    gsn_define_colormap(wks,"gsdtol")             ; choose colormap
;  gsn_define_colormap(wks,"NCV_blue_red")             ; choose colormap
;  gsn_define_colormap(wks,"MPL_gist_gray")  ; choose color map
  gsn_define_colormap(wks,"BlAqGrYeOrRevi200")  ; choose color map



   res            = True                         ; plot mods desired

   res@cnFillOn            = True                  ; color fill
   res@cnLinesOn           = True                 ; no contour lines
   res@cnInfoLabelOn       = False                ; turn off info label
   res@lbLabelBarOn        = False           ; turn off individual cb's
   res@cnLineLabelsOn      = False
   res@cnLineThicknessF      = 3.


  res@gsnScalarContour      = True              ; vectors over scalars

  res@gsnSpreadColors     = True                ; use full colormap
  res@gsnSpreadColorStart = 2
  res@gsnSpreadColorEnd   = -3


   res@cnLevelSelectionMode = "ManualLevels"  ; Manually select contours 
   res@cnMinLevelValF       = -150.     ; min
   res@cnMaxLevelValF       =  -20.     ; max
   res@cnLevelSpacingF      =   15.     ; spacing

   res@sfXCStartV = -75.      ; x stat
   res@sfXCEndV   = 75.       ; x end
   res@sfYCStartV = -75.      ; y start
   res@sfYCEndV   = 75.       ; y end

  res@sfXCStartSubsetV=-25.
  res@sfXCEndSubsetV=25.
  res@sfYCStartSubsetV=-25.
  res@sfYCEndSubsetV=25.

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
  res@vcRefMagnitudeF         = 10.0              ; make vectors larger
  res@vcRefLengthF            = 0.050            ; ref vec length
  res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.017            ; thin out vectors


  res@vfXCStartV = -75.
  res@vfXCEndV   = 75.
  res@vfYCStartV = -75.
  res@vfYCEndV   = 75.

  res@vfXCStartSubsetV=-25.
  res@vfXCEndSubsetV=25.
  res@vfYCStartSubsetV=-25.
  res@vfYCEndSubsetV=25.
     res@vcRefAnnoOn               = False

  res@vpXF              = 0.1 
  res@vpYF              = 0.9 
  res@vpWidthF          = 0.5
  res@vpHeightF         = 0.5


;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

  res@gsnMaximize  = True


   res@tiYAxisString  = "y(km)" ; y-axis title
   res@tiXAxisString  = "x(km) " ; y-axis title
;   res@tiMainString = "Fluid Elevation"  ; plot title
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 do it=0,26-1
h  = fbindirread ("shallow_h.bin",it, (/600,600/),"float")
ru  = fbindirread ("shallow_u.bin",it, (/600,600/),"float")
rv  = fbindirread ("shallow_v.bin",it, (/600,600/),"float")

   plot = gsn_csm_vector_scalar(wks,ru,rv,h,res)

 end do

end
