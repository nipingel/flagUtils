;; small script to average the power in each polarization of an ON/OFF scan pair
PRO sysFlux, sname, bm

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='Track')
;allScans = allScans[n_elements(allScans) - 9:n_elements(allScans) - 1]
allScans = allScans[4:13]
print, allScans
chan

byuBms = [1, 7, 2, 3, 4, 5, 6]

;; index allScans to select the nine scans in a 7Ptcal
scanIdx = byuBms[bm]
scanArr = [allScans[0], 1, allScans[8]]
scanArr[1] = allScans[scanIdx]

;; set channel ranges to measure statistics near HI line
a1 = 100
a2 = 145
a3 = 180
a4 = 190

;; accum buffer 0 holds power from ON_YY
;; accum buffer 1 holds power from ON_XX
;; accum buffer 2 holds power from OFF_YY
;; accum buffer 3 holds power from OFF_XX
sclear, 0
sclear, 1
sclear, 2
sclear, 3
finalMaxVal_XX = 0
finalMaxVal_YY = 0
for idx = 0, 8 DO BEGIN
    goodIdx=allScans(idx) ;this is just to avoid editing all of spencers old code
    print,'Working on scan: ', goodIdx
    get, scan=goodIdx, int=0, plnum=0; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations
    for int=0,nint-1 do begin
        for pl=0,1 do begin
           IF (idx EQ 0) || (idx EQ 8) THEN BEGIN
           	gettp, goodIdx, int=int, plnum=pl, /quiet
           	accum, pl, weight=1.0
           ENDIF 
	   if (idx GT 0) AND (idx LT 8) THEN BEGIN
		gettp, goodIdx, int=int, plnum=pl, /quiet
           	bandpass = getdata()
           	;; sum bandpass bracketing HI signal
           	currMaxVal1 = TOTAL(bandpass[a1:a2])
           	;;currMaxVal2 = TOTAL(bandpass[155:218])
           	currMaxVal2 = TOTAL(bandpass[a3:a4])
           	;; take average
           	currentMaxVal = (currMaxVal1 + currMaxVal2)/2
           	IF currentMaxVal GT finalMaxVal_YY AND pl EQ 0 THEN BEGIN
               		finalMaxVal_YY = currentMaxVal
           	ENDIF
           	IF currentMaxVal GT finalMaxVal_XX AND pl EQ 1 THEN BEGIN
               		finalMaxVal_XX = currentMaxVal
		ENDIF
	    ENDIF
        ENDFOR
    ENDFOR
ENDFOR
;ave, 0
;onPower_YYArr = getdata()
;onPower_YY1 = TOTAL(onPower_YYArr[a1:a2], /NAN)
;onPower_YY2 = TOTAL(onPower_YYArr[a3:a4], /NAN)
;onPower_YY = 0.5*(onPower_YY1 + onPower_YY2)
;copy, 0, 1
;ave, 1
;onPower_XXArr = getdata()
;onPower_XX1 = TOTAL(onPower_XXArr[a1:a2], /NAN)
;onPower_XX2 = TOTAL(onPower_XXArr[a3:a4], /NAN)
;onPower_XX = 0.5*(onPower_XX1 + onPower_XX2)
;print, onPower_XX
;copy, 0, 2
onPower_XX = finalMaxVal_XX
onPower_YY = finalMaxVal_YY
ave, 0
offPowerBandpass_YY = getdata()
copy, 0, 3
ave, 1
offPowerBandpass_XX = getdata()
copy, 0, 4
;; get powr ratio for On/Off_YY
;;subtract, 0, 2, 0
;;divide, 0, 2, 0
ave1_YY = TOTAL(offPowerBandpass_YY[a1:a2])
ave2_YY = TOTAL(offPowerBandpass_YY[a3:a4])
std1_YY = STDDEV(offPowerBandpass_YY[a1:a2])
std2_YY = STDDEV(offPowerBandpass_YY[a3:a4])
ave_YY = (ave1_YY + ave2_YY)/2
std_YY = SQRT(std1_YY^2 + std2_YY^2)/2

;; get power ratio for On/Off_XX
;;subtract, 1, 3, 0
;;divide, 0, 3, 0
ave1_XX = TOTAL(offPowerBandpass_XX[a1:a2], /NAN)
ave2_XX = TOTAL(offPowerBandpass_XX[a3:a4], /NAN)
std1_XX = STDDEV(offPowerBandpass_XX[a1:a2], /NAN)
std2_XX = STDDEV(offPowerBandpass_XX[a3:a4], /NAN)
ave_XX = (ave1_XX + ave2_XX)/2
std_XX = SQRT(std1_XX^2 + std2_XX^2)/2

offPower_YY = ave_YY
offPower_YY_Err = std_YY
offPower_XX = ave_XX
offPower_XX_Err = std_XX

print, 'Average Power in YY: ', ave_YY
print, 'Uncertainty in YY: ', std_YY
print, 'Average Power in XX: ', ave_XX
print, 'Uncertainty in XX: ', std_XX

srcFlux_XX = 22.55
srcFlux_XX_Err = 1.0017
srcFlux_YY = 22.55
srcFlux_YY_Err = 1.0017

;; calculate XX system flux
S_sys_XX = srcFlux_XX*(onPower_XX/offPower_XX - 1)^(-1)
S_sys_XX_Err = SQRT(1/(onPower_XX/offPower_XX - 1)^2*srcFlux_XX_Err^2 + ((onPower_XX*srcFlux_XX)/(onPower_XX - offPower_XX)^2)*offPower_XX_Err^2)

;; calculate YY system flux
S_sys_YY = srcFlux_YY*(onPower_YY/offPower_YY - 1)^(-1)
S_sys_YY_Err = SQRT(1/(onPower_YY/offPower_YY - 1)^2*srcFlux_YY_Err^2 + ((onPower_YY*srcFlux_YY)/(onPower_YY - offPower_YY)^2)*offPower_YY_Err^2)

print, 'XX System Flux Density: ', S_sys_XX
print, '+/-', S_sys_XX_Err
print, 'YY System Flux Density: ', S_sys_YY
print, '+/-', S_sys_YY_Err


END


