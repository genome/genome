-- Deploy result.param.software_result_id_name
-- requires: result_param

BEGIN;

CREATE INDEX result_param_id_name on result.param using btree (software_result_id, name);

COMMIT;
