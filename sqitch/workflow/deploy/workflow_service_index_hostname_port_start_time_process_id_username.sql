-- Deploy workflow_service_index_hostname_port_start_time_process_id_username
-- requires: workflow_service

BEGIN;

CREATE UNIQUE INDEX service_hostname_port_start_time_process_id_username_idx
    ON workflow.service USING btree (hostname, port, start_time, process_id, username);

COMMIT;
