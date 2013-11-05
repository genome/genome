-- Deploy result.software_result.class_name
-- requires: result_software_result

BEGIN;

CREATE INDEX sr_cname on result.software_result using btree (class_name);

COMMIT;
