-- Revert workflow_instance_index_peer_instance_id

BEGIN;

DROP INDEX workflow.instance_peer_instance_id_idx;

COMMIT;
