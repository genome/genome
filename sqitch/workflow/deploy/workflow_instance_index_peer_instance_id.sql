-- Deploy workflow_instance_index_peer_instance_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_peer_instance_id_idx ON workflow.instance USING btree (peer_instance_id);

COMMIT;
