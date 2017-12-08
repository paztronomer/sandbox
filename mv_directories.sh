for file in $(ls -d /pnfs/des/persistent/wsdiff/exp/20171119/*) ; do newname=$(echo $file | sed -r -e "s/0069/69/") ; mv $file $newname ; done
