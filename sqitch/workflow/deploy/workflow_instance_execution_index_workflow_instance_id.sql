-- Deploy workflow_instance_execution_index_workflow_instance_id
-- requires: workflow_instance_execution

BEGIN;

CREATE INDEX instance_execution_workflow_instance_id_idx
    ON workflow.instance_execution USING btree (workflow_instance_id);

COMMIT;
