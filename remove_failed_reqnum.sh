#!/bin/bash

##
## NOTE: As for July 19, 2018, the datastate.py script changes the DATA_STATE
## on PFW_ATTEMPT alone, not in ATTEMPT_STATE. The last is needed for delete
## files
##

## Variables to be used
section=db-desoper
reqnum_ls=(3493 3494)
r_table=mark2junk.txt

cnt=0
for j in ${reqnum_ls[@]}
    do
    echo "Working on reqnum=" $j
    
    ## Mark reqnum files as JUNK
    cmdlinequery.py --section $section --filename r$j_junk.aux --query "select distinct att.unitname, att.reqnum, att.attnum, att.data_state from pfw_attempt att, task t where t.id=att.task_id and t.status!=0 and att.reqnum=$j"

    cmdlinequery.py --section $section --filename r$j_null.aux --query "select distinct att.unitname, att.reqnum, att.attnum, att.data_state from pfw_attempt att, task t where t.id=att.task_id and t.status is NULL and att.reqnum=$j"
    
    cmdlinequery.py --section $section --filename r$j_good.aux --query "select distinct att.unitname, att.reqnum, att.attnum, att.data_state from pfw_attempt att, task t where t.id=att.task_id and t.status=0 and att.data_state!='JUNK' and att.reqnum=$j"
    
    ## Append to the table all but the first line
    tail -n +2 r$j_junk.aux >> $r_table   
    tail -n +2 r$j_null.aux >> $r_table
    tail -n +2 r$j_good.aux >> $r_table
    
    # Options to increase counter:
    # counter=$(expr $counter + 1)
    # ((counter++))
    ((cnt++))

    done

sleep 7
echo "Marking ALL reqnum runs as junk "
datastate.py -s $section -n JUNK -f $r_table --dbupdate

## Delete complete reqnum
for k in ${reqnum_ls[@]}
    do
    echo "Deleting reqnum=" $k
    delete_files.py --section $section --archive desar2home --reqnum $k 
    sleep 7
    done
