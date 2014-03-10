-- Deploy workflow_instance_index_parent_instance_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_parent_instance_id_idx ON workflow.instance USING btree (parent_instance_id);

COMMIT;
