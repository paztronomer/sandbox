UPDATE ATTEMPT_STATE st
    SET st.data_state='JUNK'
    WHERE st.pfw_attempt_id in (select att.id from pfw_attempt att where att.reqnum=3494);
