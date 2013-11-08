-- Deploy workflow_instance_index_name
-- requires: workflow_instance

BEGIN;

CREATE INDEX instance_name_idx ON workflow.instance USING btree (name);

COMMIT;
