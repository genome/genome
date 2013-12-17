-- Deploy workflow_instance_index_name_workflow_instance_id
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_name_workflow_instance_id_idx ON workflow.instance USING btree (name, workflow_instance_id);

COMMIT;
