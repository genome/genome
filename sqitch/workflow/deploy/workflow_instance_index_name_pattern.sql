-- Deploy workflow_instance_index_name_pattern
-- requires: workflow_instance

BEGIN;

CREATE INDEX workflow_instance_name_varchar_pattern_ops_index
    ON workflow.instance USING btree (name varchar_pattern_ops);

COMMIT;
