PRO PAF_combineDither, sname, fileout=outFile, bm=bm, sess = sess

fileout, outfile

;; read in session with lower LO
filein, '/home/scratch/npingel/FLAG/data/AGBT17B_360/AGBT17B_360_03/NGC4258_Field/rawBeamFormedData/AGBT17B_360_03_NGC4258_Field_Beam' + bm + '_edge_chanBlank_freqUpdate_Baselined.fits'
;;goodScans_1 is for 17B_360_03 (lower LO setting)
goodScans_1 = get_scan_numbers(source=sname,procedure='*Map')
goodScans_1 = goodScans[0:31]


;; read in session with higher LO
filein, 'NGC6946_edge_13_Beam' + bm + '_channelShift.fits'
;; select all scans of a source that are mapping scans
;; goodScans_2 is for higher LO setting data
goodScans_2 = get_scan_numbers(source=sname,procedure='*Map')

freeze
;;chan
cnt = 0
for idx = 0,N_ELEMENTS(goodScans_2)-1 DO BEGIN
;;for idx = 0,1 DO BEGIN
    sclear, 0
    sclear, 1
    goodScan_1 = goodScans_1(cnt)
    goodScan_2 = goodScans_2(idx)
    print,'Working on scan: ', goodScan_2
    get, scan=goodScan_2, int=0, plnum=1; call scan to get the number of intergrations
    ;; get scan info
    info=scan_info(goodScan_2) ;; so to get 63 integrations as oppposed to 70 in session 12)
    nint=info.n_integrations
    for int=0,nint-1 do begin 
        for pl=0,1 do begin
           get, scan=goodScan_2, int=int, plnum=pl
           dcshift, !g.s[0], 16 ;; shift higher LO data
           accum, pl
           filein, 'AGBT17B_360_0' + sess + '_NGC4258_Field_Beam' + bm + '_chanBlank_edge_Baselined.fits'
           get, scan=goodScan_1, int=int, plnum=pl
           accum, pl
           ave, pl
           keep
           filein, 'AGBT17B_360_' + sname + '_Comb_03_0' + sess + '_edge_Beam' + bm + '.fits'
        ENDFOR
    ENDFOR
    cnt = cnt + 1
    IF cnt EQ 32 THEN BEGIN
    	cnt = 0
    ENDIF
ENDFOR
unfreeze
END


