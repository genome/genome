-- Deploy result.param.software_result_id
-- requires: result_param

BEGIN;

CREATE INDEX srp_sri2 on result.param using btree (software_result_id);

COMMIT;
