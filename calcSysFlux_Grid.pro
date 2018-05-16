;; small script to average the power in each polarization of an ON/OFF scan pair
PRO calcSysFlux_Grid, sname

freeze
chan
;; select all scans of a source that are mapping scans
gridScans=get_scan_numbers(source=sname,procedure='*Map')

a0 = 125
a1 = 134
a2 = 155
a3 = 180

;; get all scans from the Track OFF position
offScans=get_scan_numbers(source=sname, procedure='Track')
print, offScans

;; find max value of mapping scans
finalMaxVal_XX = 0
finalMaxVal_YY = 0
scanNum_XX = 0
intNum_XX = 0
scanNum_YY = 0
intNum_YY = 0

srcFlux_XX = 22.2294 
srcFlux_XX_Err = 1.0017
srcFlux_YY = 22.2294
srcFlux_YY_Err = 1.0017
for idx = 0,N_ELEMENTS(gridScans)-1 DO BEGIN
    goodIdx=gridScans[idx] ;this is just to avoid editing all of spencers old code
    print,'Working on scan: ', goodIdx
    get, scan=goodIdx, int=0, plnum=0; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations
    for int=0,nint-1 do begin
        for pl=0,1 do begin
           gettp, goodIdx, int=int, plnum=pl, /quiet
           bandpass = getdata()
           ;; sum bandpass bracketing HI signal
           currMaxVal1 = TOTAL(bandpass[a0:a1])
	   currMaxVal2 = TOTAL(bandpass[a2:a3])
           ;; take average
           currentMaxVal = (currMaxVal1 + currMaxVal2)/2
           IF currentMaxVal GT finalMaxVal_YY AND pl EQ 0 THEN BEGIN
               finalMaxVal_YY = currentMaxVal
               scanNum_YY = goodIdx
               intNum_YY = int
           ENDIF
           IF currentMaxVal GT finalMaxVal_XX AND pl EQ 1 THEN BEGIN
               finalMaxVal_XX = currentMaxVal
               scanNum_XX = goodIdx
               intNum_XX = int          
           ENDIF
           ENDFOR
    ENDFOR
ENDFOR

print, 'Max Total Power (XX): ', finalMaxVal_XX
print, 'Scan (XX): ', scanNum_XX
print, 'Int (XX): ', intNum_XX

print, 'Max Total Power (YY): ', finalMaxVal_YY
print, 'Scan (YY): ', scanNum_YY
print, 'Int (YY): ', intNum_YY


onPower_XX = finalMaxVal_XX
onPower_YY = finalMaxVal_YY
finalMaxValOffYY = 0
finalMaxValOffXX = 0
for idx = 0, N_ELEMENTS(offScans) -1 DO BEGIN 
    goodIdx=offScans[idx] ;this is just to avoid editing all of spencers old code
    print,'Working on scan: ', goodIdx
    get, scan=goodIdx, int=0, plnum=0; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations
    for int=0,nint-1 do begin
        for pl=0,1 do begin
        	gettp, goodIdx, int=int, plnum=pl, /quiet
		bandpass = getdata()
           	;; sum bandpass bracketing HI signal
           	currMaxVal1 = TOTAL(bandpass[a0:a1])
           	currMaxVal2 = TOTAL(bandpass[a2:a3])
           	;; take average
           	currentMaxVal = (currMaxVal1 + currMaxVal2)/2
           	IF finalMaxValOffYY LT finalMaxValOffYY AND pl EQ 0 THEN BEGIN
               		finalMaxValYY = currentMaxVal
               		scanNum_YY = goodIdx
               		intNum_YY = int
			copy, 0, 1
           	ENDIF
           	IF currentMaxVal LT finalMaxVal_XX AND pl EQ 1 THEN BEGIN
               		finalMaxValXX = currentMaxVal
               		scanNum_XX = goodIdx
               		intNum_XX = int
			copy, 0, 2
           ENDIF
            ;;accum, pl
        ENDFOR
    ENDFOR
ENDFOR

;;get mean YY Off bandpass
;;ave, 0
copy, 1, 0
offBandpass = getdata()
offPower_YY1 = TOTAL(offBandpass[a0:a1])
offPower_YY2 = TOTAL(offBandpass[a2:a3])
offPower_YY_Err1 = STDDEV(offBandpass[a0:a1])
offPower_YY_Err2 = STDDEV(offBandpass[a2:a3])
offPower_YY = (offPower_YY1 + offPower_YY2)/2
offPower_YY_Err = SQRT(offPower_YY_Err1^2 + offPower_YY_Err2^2)/2

;;get mean XX Off bandpass
;;ave, 1
copy, 2, 0
offBandpass = getdata()
offPower_XX1 = TOTAL(offBandpass[a0:a1])
offPower_XX2 = TOTAL(offBandpass[a2:a3])
offPower_XX_Err1 = STDDEV(offBandpass[a0:a1])
offPower_XX_Err2 = STDDEV(offBandpass[a2:a3])
offPower_XX = (offPower_XX1 + offPower_XX2)/2
offPower_XX_Err = SQRT(offPower_XX_Err1^2 + offPower_XX_Err2^2)/2
;; calculate XX system flux
S_sys_XX = srcFlux_XX*(onPower_XX/offPower_XX - 1)^(-1)
S_sys_XX_Err = SQRT(1/(onPower_XX/offPower_XX - 1)^2*srcFlux_XX_Err^2 + ((onPower_XX*srcFlux_XX)/(onPower_XX - offPower_XX)^2)*offPower_XX_Err^2)

;; calculate YY system flux
S_sys_YY = srcFlux_YY*(onPower_YY/offPower_YY - 1)^(-1)
S_sys_YY_Err = SQRT(1/(onPower_YY/offPower_YY - 1)^2*srcFlux_YY_Err^2 + ((onPower_YY*srcFlux_YY)/(onPower_YY - offPower_YY)^2)*offPower_YY_Err^2)

print, 'On Power: ', onPower_XX
print, 'Mean Off Power (XX)', offPower_XX
print, '+.-',offPower_XX_Err
print, 'XX System Flux Density: ', S_sys_XX
print, '+/-', S_sys_XX_Err

print, 'On Power: ', onPower_YY
print, 'Mean Off Power (YY)', offPower_YY
print, '+.-',offPower_YY_Err
print, 'YY System Flux Density: ', S_sys_YY
print, '+/-', S_sys_YY_Err

END


