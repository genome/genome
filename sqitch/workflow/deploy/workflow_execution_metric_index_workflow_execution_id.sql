-- Deploy workflow_execution_metric_index_workflow_execution_id
-- requires: workflow_execution_metric

BEGIN;

CREATE INDEX execution_metric_workflow_execution_id_idx ON workflow.execution_metric
        USING btree (workflow_execution_id);

COMMIT;
