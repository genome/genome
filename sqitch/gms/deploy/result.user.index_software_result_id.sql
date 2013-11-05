-- Deploy result.user.software_result_id
-- requires: result_user

BEGIN;

CREATE INDEX sru_rid_i on result."user" using btree (software_result_id);

COMMIT;
