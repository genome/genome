-- Deploy workflow_instance_execution_index_workflow_execution_id
-- requires: workflow_instance_execution

BEGIN;

CREATE UNIQUE INDEX instance_execution_workflow_execution_id_idx ON workflow.instance_execution
        USING btree (workflow_execution_id);

COMMIT;
