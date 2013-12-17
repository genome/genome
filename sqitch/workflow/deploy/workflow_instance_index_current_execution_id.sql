-- Deploy workflow_instance_index_current_execution_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_current_execution_id_idx ON workflow.instance USING btree (current_execution_id);

COMMIT;
