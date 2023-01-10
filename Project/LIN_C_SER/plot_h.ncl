;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
;************************************************
begin
;  wks = gsn_open_wks("ncgm" ,"shallow")             ; open a ps file
  wks = gsn_open_wks("pdf" ,"shallow")             ; open a ps file
;  gsn_define_colormap(wks,"gsltod")             ; choose colormap
;    gsn_define_colormap(wks,"gsdtol")             ; choose colormap
;  gsn_define_colormap(wks,"NCV_blue_red")             ; choose colormap
  gsn_define_colormap(wks,"MPL_gist_gray")  ; choose color map


   res            = True                         ; plot mods desired

   res@cnFillOn            = True                  ; color fill
   res@cnLinesOn           = False                 ; no contour lines
   res@cnInfoLabelOn       = False                ; turn off info label
   res@lbLabelBarOn        = False           ; turn off individual cb's
   res@cnLineLabelsOn      = False

   res@cnLevelSelectionMode = "ManualLevels"  ; Manually select contours 
   res@cnMinLevelValF       = -0.9     ; min
   res@cnMaxLevelValF       =  0.9     ; max
   res@cnLevelSpacingF      =  0.2     ; spacing

   res@sfXCStartV = -30.      ; x stat
   res@sfXCEndV   = 30.       ; x end
   res@sfYCStartV = -30.      ; y start
   res@sfYCEndV   = 30.       ; y end

   res@vpXF              = 0.2  ; left viewport 
   res@vpYF              = 0.8  ; top viewport 
   res@vpWidthF          = 0.6
   res@vpHeightF         = 0.6

   res@tiYAxisString  = "y(km)" ; y-axis title
   res@tiXAxisString  = "x(km) " ; y-axis title
   res@tiMainString = "Fluid Elevation"  ; plot title
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; do it=9,9
  do it=0,9
   h  = fbindirread ("shallow.bin",it, (/600,600/),"float")  ; read binary record it
   plot = gsn_csm_contour(wks,h,res)         ; plot
 end do

end
