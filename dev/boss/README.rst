=================
specter/dev/boss/
=================

This directory is for *code* related to using specter with BOSS.
It is *not* for data files such as PSFs -- keep this directory lean
and light so that it doesn't bog down non-BOSS users of specter.

::

    cd $SPECTER_DIR/data

    specter -i calib/arc/arc_lines.fits -o img-boss-arc.fits \
        -I dark/boss/dark-900-r1-00124292.fits \
        -r 0,20 \
        -p boss/pixpsf-r1-00140299.fits \
        -t bigboss/designs/20120827difdet/bbthru-R.fits

Hmm.  Need BOSS throughput model to get this right...
