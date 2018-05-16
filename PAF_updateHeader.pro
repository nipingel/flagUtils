PRO PAF_updateHeader, sname, fileout=outFile

;; select all scans of a source that are mapping scans
allScans=get_scan_numbers(source=sname,procedure='*Map')
fileout, outfile

;; remove badscans
IF NOT KEYWORD_SET(badScans) THEN badScans=[0]
goodScans=setdifference(allScans,badScans)

freeze
chan
for idx = 0,N_ELEMENTS(goodScans)-1 DO BEGIN
    goodIdx=goodScans(idx) ;this is just to avoid editing all of spencers old code

    print,'Working on scan: ', goodIdx
    get, scan=goodIdx, int=0, plnum=0; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodIdx)
    nint=info.n_integrations
    for int=0,nint-1 do begin 
        for pl=0,1 do begin
           get, scan=goodIdx, int=int, plnum=pl
           !g.s[0].reference_frequency = !g.s[0].reference_frequency - 151.59*1e3
           !g.s[0].center_frequency = !g.s[0].center_frequency - 151.59*1e3
           data = getdata()
           for i=0,98 do begin
               strtChan = 28+i*32
               endChan = strtChan + 6
               data[strtChan:endChan] = !Values.F_NAN
           ENDFOR
           setdata, data
           keep
        ENDFOR
    ENDFOR
ENDFOR
unfreeze
END


