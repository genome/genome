-- Deploy workflow_instance_index_parent_execution_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_parent_execution_id_idx ON workflow.instance USING btree (parent_execution_id);

COMMIT;
