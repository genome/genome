-- Verify workflow_service_index_hostname_port_start_time_process_id_username

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'service_hostname_port_start_time_process_id_username_idx';

ROLLBACK;
