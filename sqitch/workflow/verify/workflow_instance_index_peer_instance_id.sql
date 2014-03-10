-- Verify workflow_instance_index_peer_instance_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instance_peer_instance_id_idx';

ROLLBACK;
