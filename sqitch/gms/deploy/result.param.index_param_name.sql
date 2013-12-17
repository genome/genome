-- Deploy result.param.param_name
-- requires: result_param

BEGIN;

CREATE INDEX srp_pn2 on result.param using btree (param_name);

COMMIT;
