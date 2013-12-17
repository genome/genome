-- Deploy workflow_instance_index_workflow_plan_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_workflow_plan_id_idx ON workflow.instance USING btree (workflow_plan_id);

COMMIT;
