for f in *xtalked.fits; do mv "$f" "`echo $f | sed s/initial_string_pattern/final_string_pattern/`"; done
