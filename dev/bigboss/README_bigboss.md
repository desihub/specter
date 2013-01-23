specter/dev/bigboss/
====================

This directory is for *code* related to using specter with BigBOSS.
It is *not* for data files such as PSFs -- keep this directory lean
and light so that it doesn't bog down non-BigBOSS users of specter.

### Example commands for quick reference ###

<pre>
cd $SPECTER_DIR/data/
specter -i sky/sky-uves.fits -o blat.fits -r 0,20 -w 7900,8000 --extra \
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
</pre>
