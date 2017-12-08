for file in $(ls /pnfs/des/persistent/wsdiff/exp/20171119/*/*.fz) ; do newname=$(echo $file | sed -r -e "s/_c([0-9])/_\1/" -e "s/immasked/immask/") ; mv $file $newname ; done
