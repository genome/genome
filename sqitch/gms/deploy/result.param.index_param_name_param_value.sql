-- Deploy result.param.param_name_param_value
-- requires: result_param

BEGIN;

CREATE INDEX srp_pnpv2 on result.param using btree (param_name, param_value);

COMMIT;
