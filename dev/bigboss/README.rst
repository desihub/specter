====================
specter/dev/bigboss/
====================

Introduction
------------

This directory is for *code* related to using specter with BigBOSS.
It is *not* for data files such as PSFs -- keep this directory lean
and light so that it doesn't bog down non-BigBOSS users of specter.

Example commands for quick reference
------------------------------------

Examples::

    #- Arc
    cd $SPECTER_DIR/data/
    specter -i calib/arc/arc_continuum.fits -o blat.fits \
        -r 0,20 -w 7900,8000 --extra \
        -p bigboss/designs/20120827difdet/bbpsf-I.fits

    specter -i calib/arc/arc_lines.fits -o img-arc.fits \
        -I dark/boss/dark-900-r1-00124292.fits \
        -r 0,20 --noise \
        -p bigboss/designs/20120827difdet/bbpsf-I.fits

    #- BBspecsim ELG grid
    cd $SPECTER_DIR/data/
    specter -i bigboss/bbspecsim/elg_grid.fits.gz -o blat.fits \
        -r 300,325 -w 8200,8575 \
        -p bigboss/designs/20120827difdet/bbpsf-I.fits

    #- Sky
    cd $SPECTER_DIR/data/
    specter -i sky/sky-uves.fits -o img-sky.fits \
        -r 0,49 -w 8300,8700 \
        --noise \
        -I dark/boss/dark-900-r1-00124292.fits \
        -p bigboss/designs/20120827difdet/bbpsf-I.fits

    cd $SPECTER_DIR/data/bigboss/designs/20120827difdet/
    $SPECTER_DIR/dev/bigboss/spots2psf.py \
        -i /data/bigboss/sim/spots/BB_SPEC_20120827difdet/Blue/ \
        -t bbthru-B.fits -o bbpsf-B.fits

    $SPECTER_DIR/dev/bigboss/spots2psf.py \
        -i /data/bigboss/sim/spots/BB_SPEC_20120827difdet/Red_250/ \
        -t bbthru-R.fits -o bbpsf-R.fits

    $SPECTER_DIR/dev/bigboss/spots2psf.py \
        -i /data/bigboss/sim/spots/BB_SPEC_20120827difdet/NIR_500/ \
        -t bbthru-I.fits -o bbpsf-I.fits

    cd $SPECTER_DIR/data/
    specter -i sky/sky-uves.fits -o blat.fits -r 0,1 \
            -p bigboss/designs/20120827difdet/bbpsf-I-500um.fits \
            -t bigboss/designs/20120827difdet/bbthru-I.fits

    cd $SPECTER_DIR/data/
    specter -i sky/sky-uves.fits -o blat.fits -r 0,1 -w 7900,8000 \
            -p bigboss/designs/20120827difdet/blat.fits

    cd $SPECTER_DIR/dev/bigboss/
    ./spots2psf.py -i $BBSPECSIM_DIR/designs/BB_SPEC_20120428difdet/spots/Blue/ \
        -t bbthru-B.fits -o blat.fits

    ./spots2psf.py -i /data/bigboss/sim/spots/BB_SPEC_20120827difdet/Blue/ \
        -t bbthru-B.fits -o blat.fits
