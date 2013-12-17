-- Deploy result.param.name
-- requires: result_param

BEGIN;

CREATE INDEX result_param_name_index on result.param using btree (name);

COMMIT;
